/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/**
 * Copyright (c) 2011-2015  Tsinghua University, P.R.China.
 *
 * This file is part of ndnSIM. See AUTHORS for complete list of ndnSIM authors and
 * contributors.
 *
 * ndnSIM is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * ndnSIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * ndnSIM, e.g., in COPYING.md file.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @author Xiaoke Jiang <shock.jiang@gmail.com>
 **/

#include "ndn-consumer-versioning-geometric.hpp"

#include <math.h>

NS_LOG_COMPONENT_DEFINE("ndn.ConsumerVersioningGeometric");

namespace ns3 {
namespace ndn {

NS_OBJECT_ENSURE_REGISTERED(ConsumerVersioningGeometric);

TypeId
ConsumerVersioningGeometric::GetTypeId(void)
{
  static TypeId tid =
    TypeId("ns3::ndn::ConsumerVersioningGeometric")
      .SetGroupName("Ndn")
      .SetParent<ConsumerCbr>()
      .AddConstructor<ConsumerVersioningGeometric>()

      .AddAttribute("NumberOfVersions", "Number of the Versions in total", StringValue("10"),
                    MakeUintegerAccessor(&ConsumerVersioningGeometric::SetNumberOfVersions,
                                         &ConsumerVersioningGeometric::GetNumberOfVersions),
                    MakeUintegerChecker<uint32_t>())

      .AddAttribute("p", "parameter of power", StringValue("0.7"),
                    MakeDoubleAccessor(&ConsumerVersioningGeometric::SetP,
                                       &ConsumerVersioningGeometric::GetP),
                    MakeDoubleChecker<double>());

  return tid;
}

ConsumerVersioningGeometric::ConsumerVersioningGeometric()
  : m_N(10)
  , m_p(0.7)
  , m_seqRng(CreateObject<UniformRandomVariable>())
{
  // SetNumberOfVersions is called by NS-3 object system during the initialization
}

ConsumerVersioningGeometric::~ConsumerVersioningGeometric()
{
}

void
ConsumerVersioningGeometric::SetNumberOfVersions(uint32_t numOfVersions)
{
  m_N = numOfVersions;

  NS_LOG_DEBUG(m_p << " and " << m_N);

  m_Pcum = std::vector<double>(m_N + 1);

  m_Pcum[0] = 0.0;
  for (uint32_t i = 1; i <= m_N; i++) {
    m_Pcum[i] = m_Pcum[i - 1] + m_p * std::pow(1.0-m_p, i-1);
  }

  for (uint32_t i = 1; i <= m_N; i++) {
    m_Pcum[i] = m_Pcum[i] / m_Pcum[m_N];
    NS_LOG_LOGIC("Cumulative probability [" << i << "]=" << m_Pcum[i]);
  }
}

uint32_t
ConsumerVersioningGeometric::GetNumberOfVersions() const
{
  return m_N;
}

void
ConsumerVersioningGeometric::SetP(double p)
{
  m_p = p;
  SetNumberOfVersions(m_N);
}

double
ConsumerVersioningGeometric::GetP() const
{
  return m_p;
}

void
ConsumerVersioningGeometric::SendPacket()
{
  if (!m_active)
    return;

  NS_LOG_FUNCTION_NOARGS();

  uint32_t seq = std::numeric_limits<uint32_t>::max(); // invalid

  // std::cout << Simulator::Now ().ToDouble (Time::S) << "s max -> " << m_seqMax << "\n";

  while (m_retxSeqs.size()) {
    seq = *m_retxSeqs.begin();
    m_retxSeqs.erase(m_retxSeqs.begin());

    // NS_ASSERT (m_seqLifetimes.find (seq) != m_seqLifetimes.end ());
    // if (m_seqLifetimes.find (seq)->time <= Simulator::Now ())
    //   {

    //     NS_LOG_DEBUG ("Expire " << seq);
    //     m_seqLifetimes.erase (seq); // lifetime expired. Trying to find another unexpired
    //     sequence number
    //     continue;
    //   }
    NS_LOG_DEBUG("=interest seq " << seq << " from m_retxSeqs");
    break;
  }

  if (seq == std::numeric_limits<uint32_t>::max()) // no retransmission
  {
    if (m_seqMax != std::numeric_limits<uint32_t>::max()) {
      if (m_seq >= m_seqMax) {
        return; // we are totally done
      }
    }

    seq = ConsumerVersioningGeometric::GetNextSeq();
    m_seq++;
  }

  // std::cout << Simulator::Now ().ToDouble (Time::S) << "s -> " << seq << "\n";

  //
  shared_ptr<Name> nameWithVersion = make_shared<Name>(m_interestName);
#if 0  /* SEGV occurs without some sequence number? */
  nameWithVersion->appendVersion(seq);
#else
  nameWithVersion->appendSequenceNumber(seq);
#endif
  //

  shared_ptr<Interest> interest = make_shared<Interest>();
  interest->setNonce(m_rand->GetValue(0, std::numeric_limits<uint32_t>::max()));
  interest->setName(*nameWithVersion);

  // NS_LOG_INFO ("Requesting Interest: \n" << *interest);
  NS_LOG_INFO("> Interest for " << seq << ", Total: " << m_seq << ", face: " << m_face->getId());
  NS_LOG_DEBUG("Trying to add " << seq << " with " << Simulator::Now() << ". already "
                                << m_seqTimeouts.size() << " items");

  m_seqTimeouts.insert(SeqTimeout(seq, Simulator::Now()));
  m_seqFullDelay.insert(SeqTimeout(seq, Simulator::Now()));

  m_seqLastDelay.erase(seq);
  m_seqLastDelay.insert(SeqTimeout(seq, Simulator::Now()));

  m_seqRetxCounts[seq]++;

  m_rtt->SentSeq(SequenceNumber32(seq), 1);

  m_transmittedInterests(interest, this, m_face);
  m_appLink->onReceiveInterest(*interest);

  ConsumerVersioningGeometric::ScheduleNextPacket();
}

uint32_t
ConsumerVersioningGeometric::GetNextSeq()
{
  uint32_t content_index = 1; //[1, m_N]
  double p_sum = 0;

  double p_random = m_seqRng->GetValue();
  while (p_random == 0) {
    p_random = m_seqRng->GetValue();
  }
  // if (p_random == 0)
  NS_LOG_LOGIC("p_random=" << p_random);
  for (uint32_t i = 1; i <= m_N; i++) {
    p_sum = m_Pcum[i]; // m_Pcum[i] = m_Pcum[i-1] + p[i], p[0] = 0;   e.g.: p_cum[1] = p[1],
                       // p_cum[2] = p[1] + p[2]
    if (p_random <= p_sum) {
      content_index = i;
      break;
    } // if
  }   // for
  // content_index = 1;
  NS_LOG_DEBUG("RandomNumber=" << content_index);
  return content_index;
}

void
ConsumerVersioningGeometric::ScheduleNextPacket()
{

  if (m_firstTime) {
    m_sendEvent = Simulator::Schedule(Seconds(0.0), &ConsumerVersioningGeometric::SendPacket, this);
    m_firstTime = false;
  }
  else if (!m_sendEvent.IsRunning())
    m_sendEvent = Simulator::Schedule((m_random == 0) ? Seconds(1.0 / m_frequency)
                                                      : Seconds(m_random->GetValue()),
                                      &ConsumerVersioningGeometric::SendPacket, this);
}

} /* namespace ndn */
} /* namespace ns3 */
