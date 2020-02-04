/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License version 2 as
* published by the Free Software Foundation;
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program; if not, write to the Free Software
* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#undef PGO_TRAINING
#define PATH_TO_PGO_CONFIG "path_to_pgo_config"

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <random>
#include <time.h> 
#include "ns3/core-module.h"
#include "ns3/qbb-helper.h"
#include "ns3/point-to-point-helper.h"
#include "ns3/applications-module.h"
#include "ns3/internet-module.h"
#include "ns3/global-route-manager.h"
#include "ns3/ipv4-static-routing-helper.h"
//#include "ns3/broadcom-node.h"
#include "ns3/packet.h"
#include "ns3/error-model.h"
#include <ns3/rdma.h>
#include <ns3/rdma-client.h>
#include <ns3/rdma-client-helper.h>
#include <ns3/rdma-driver.h>
#include <ns3/switch-node.h>
#include <ns3/sim-setting.h>
#include <vector>
#include <set>
#include <unordered_map>

using namespace ns3;
using namespace std;

NS_LOG_COMPONENT_DEFINE("GENERIC_SIMULATION");

uint32_t cc_mode = 1;
bool enable_qcn = true, use_dynamic_pfc_threshold = true;
uint32_t packet_payload_size = 1000, l2_chunk_size = 0, l2_ack_interval = 0;
double pause_time = 5, simulator_stop_time = 3.01;
std::string data_rate, link_delay, topology_file, flow_file, trace_file, trace_output_file;
std::string fct_output_file = "fct.txt";
std::string pfc_output_file = "pfc.txt";



double app_start_time = 1.0;

double alpha_resume_interval = 55, rp_timer, ewma_gain = 1 / 16;
double rate_decrease_interval = 4;
uint32_t fast_recovery_times = 5;
std::string rate_ai, rate_hai, min_rate = "100Mb/s";
std::string dctcp_rate_ai = "1000Mb/s";

bool clamp_target_rate = false, l2_back_to_zero = false;
double error_rate_per_link = 0.0;
uint32_t has_win = 1;
uint32_t global_t = 1;
uint32_t mi_thresh = 5;
bool var_win = false, fast_react = true;
bool multi_rate = true;
bool sample_feedback = false;
double u_target = 0.95;
uint32_t int_multi = 1;
bool rate_bound = true;

uint32_t ack_high_prio = 0;
uint64_t link_down_time = 0;
uint32_t link_down_A = 0, link_down_B = 0;

uint32_t enable_trace = 1;

uint32_t buffer_size = 16;

uint32_t qCnt = 8;

uint32_t qlen_dump_interval = 100000000, qlen_mon_interval = 100;
uint64_t qlen_mon_start = 2000000000, qlen_mon_end = 2100000000;
string qlen_mon_file;

double sigma = 2.0;

uint32_t datarate = 100;

double app_stop_time = 2.05;

unordered_map<uint64_t, uint32_t> rate2kmax, rate2kmin;
unordered_map<uint64_t, double> rate2pmax;

bool w_dctcp = false;
bool w_3 = false;
bool conga = false;

int min_size = 100;

std::vector<uint64_t> enterprise_size_conga{
    250, 270, 286, 297, 312, 327, 342, 356, 368, 383, 404, 427, 447,
    469, 490, 512, 536, 561, 587, 612, 636, 651, 674, 696, 708, 729,
    753, 782, 820, 857, 900, 945, 990, 1036, 1086, 1140, 1200, 1260,
    1311, 1353, 1393, 1452, 1485, 1510, 1540, 1597, 1658, 1721, 1784,
    1873, 1978, 2079, 2160, 2267, 2395, 2536, 2688, 2862, 2960, 3040,
    3163, 3390, 3603, 3836, 4124, 4406, 4540, 4673, 4952, 5260, 5625,
    5938, 6088, 6506, 6976, 7430, 7800, 8332, 8849, 9477, 10344, 11256,
    12396, 13766, 15370, 17417, 19884, 22822, 26285, 30973, 37152, 45274,
    55796, 70997, 94050, 131474, 189792, 344538, 898190, 2221957, 2624209,
    3114668, 3849851, 4904979, 6438497, 8608848, 13392072, 23110447,
    40607470, 44466961, 48564256, 53520618, 59661017, 67292218, 76720988,
    92984411, 121274824, 168707183, 176569410, 184086600, 195504461,
    211698244, 230332119, 249070260, 293849828, 361127656, 440186178,
    442950536, 449947224, 466035021, 482515012, 510665738, 561765744,
    669621992, 738727896, 1022256099, 1773379248};

std::vector<double> enterprise_prob_conga{
    0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1,
    0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21,
    0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32,
    0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43,
    0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54,
    0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65,
    0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76,
    0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87,
    0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98,
    0.99, 0.991, 0.992, 0.993, 0.994, 0.995, 0.996, 0.997, 0.998,
    0.999, 0.9991, 0.9992, 0.9993, 0.9994, 0.9995, 0.9996, 0.9997,
    0.9998, 0.9999, 0.99991, 0.99992, 0.99993, 0.99994, 0.99995,
    0.99996, 0.99997, 0.99998, 0.99999, 0.999991, 0.999992, 0.999993,
    0.999994, 0.999995, 0.999996, 0.999997, 0.999998, 0.999999, 1.0};

//DCTCP
std::vector<uint64_t> enterprise_size_dctcp{0, 100, 10000, 20000, 30000, 50000, 80000, 200000, 1000000, 2000000, 5000000, 10000000, 30000000};
std::vector<double> enterprise_prob_dctcp{0.0, 0.0, 0.15, 0.2, 0.3, 0.4, 0.53, 0.6, 0.7, 0.8, 0.9, 0.97, 1.0};


std::vector<uint64_t> enterprise_size;
std::vector<double> enterprise_prob;

std::set<double> probs;
std::unordered_map<double, int> size_map;

bool poison_distr = true;
double utilization_factor = 0.6;

bool create_incast = false;
uint32_t packet_size_incast = 20;
 uint32_t num_sources = 200;

bool inter_dc_traffic = true;
uint64_t nic_rate;

uint64_t maxRtt, maxBdp;

struct Interface{
	uint32_t idx;
	bool up;
	uint64_t delay;
	uint64_t bw;

	Interface() : idx(0), up(false){}
};
map<Ptr<Node>, map<Ptr<Node>, Interface> > nbr2if;
// Mapping destination to next hop for each node: <node, <dest, <nexthop0, ...> > >
map<Ptr<Node>, map<Ptr<Node>, vector<Ptr<Node> > > > nextHop;
map<Ptr<Node>, map<Ptr<Node>, uint64_t> > pairDelay;
map<Ptr<Node>, map<Ptr<Node>, uint64_t> > pairTxDelay;
map<Ptr<Node>, map<Ptr<Node>, uint64_t> > pairBw;
map<Ptr<Node>, map<Ptr<Node>, uint64_t> > pairBdp;
map<Ptr<Node>, map<Ptr<Node>, uint64_t> > pairRtt;

void qp_finish(FILE* fout, Ptr<RdmaQueuePair> q){
	//fprintf(fout, "%lu QP complete\n", Simulator::Now().GetTimeStep());
	fprintf(fout, "%08x %08x %u %u %lu %lu %lu\n", q->sip.Get(), q->dip.Get(), q->sport, q->dport, q->m_size, q->startTime.GetTimeStep(), (Simulator::Now() - q->startTime).GetTimeStep());
	fflush(fout);
}

void get_pfc(FILE* fout, Ptr<QbbNetDevice> dev, uint32_t type){
	fprintf(fout, "%lu %u %u %u %u\n", Simulator::Now().GetTimeStep(), dev->GetNode()->GetId(), dev->GetNode()->GetNodeType(), dev->GetIfIndex(), type);
}

struct QlenDistribution{
	vector<uint32_t> cnt; // cnt[i] is the number of times that the queue len is i KB

	void add(uint32_t qlen){
		uint32_t kb = qlen / 1000;
		if (cnt.size() < kb+1)
			cnt.resize(kb+1);
		cnt[kb]++;
	}
};
map<uint32_t, map<uint32_t, QlenDistribution> > queue_result;
void monitor_buffer(FILE* qlen_output, NodeContainer *n){
	for (uint32_t i = 0; i < n->GetN(); i++){
		if (n->Get(i)->GetNodeType() == 1){ // is switch
			Ptr<SwitchNode> sw = DynamicCast<SwitchNode>(n->Get(i));
			if (queue_result.find(i) == queue_result.end())
				queue_result[i];
			for (uint32_t j = 1; j < sw->GetNDevices(); j++){
				uint32_t size = 0;
				for (uint32_t k = 0; k < qCnt; k++)
					size += sw->m_mmu->egress_bytes[j][k];
				queue_result[i][j].add(size);
			}
		}
	}
	if (Simulator::Now().GetTimeStep() % qlen_dump_interval == 0){
		fprintf(qlen_output, "time: %lu\n", Simulator::Now().GetTimeStep());
		for (auto &it0 : queue_result)
			for (auto &it1 : it0.second){
				fprintf(qlen_output, "%u %u", it0.first, it1.first);
				auto &dist = it1.second.cnt;
				for (uint32_t i = 0; i < dist.size(); i++)
					fprintf(qlen_output, " %u", dist[i]);
				fprintf(qlen_output, "\n");
			}
		fflush(qlen_output);
	}
	if (Simulator::Now().GetTimeStep() < qlen_mon_end)
		Simulator::Schedule(NanoSeconds(qlen_mon_interval), &monitor_buffer, qlen_output, n);
}

void CalculateRoute(Ptr<Node> host){
	// queue for the BFS.
	vector<Ptr<Node> > q;
	// Distance from the host to each node.
	map<Ptr<Node>, int> dis;
	map<Ptr<Node>, uint64_t> delay;
	map<Ptr<Node>, uint64_t> txDelay;
	map<Ptr<Node>, uint64_t> bw;
	// init BFS.
	q.push_back(host);
	dis[host] = 0;
	delay[host] = 0;
	txDelay[host] = 0;
	bw[host] = 0xfffffffffffffffflu;
	// BFS.
	for (int i = 0; i < (int)q.size(); i++){
		Ptr<Node> now = q[i];
		int d = dis[now];
		for (auto it = nbr2if[now].begin(); it != nbr2if[now].end(); it++){
			// skip down link
			if (!it->second.up)
				continue;
			Ptr<Node> next = it->first;
			// If 'next' have not been visited.
			if (dis.find(next) == dis.end()){
				dis[next] = d + 1;
				delay[next] = delay[now] + it->second.delay;
				txDelay[next] = txDelay[now] + packet_payload_size * 1000000000lu * 8 / it->second.bw;
				bw[next] = std::min(bw[now], it->second.bw);
				// we only enqueue switch, because we do not want packets to go through host as middle point
				if (next->GetNodeType() == 1)
					q.push_back(next);
			}
			// if 'now' is on the shortest path from 'next' to 'host'.
			if (d + 1 == dis[next]){
				nextHop[next][host].push_back(now);
			}
		}
	}
	for (auto it : delay)
		pairDelay[it.first][host] = it.second;
	for (auto it : txDelay)
		pairTxDelay[it.first][host] = it.second;
	for (auto it : bw)
		pairBw[it.first][host] = it.second;
}

void CalculateRoutes(NodeContainer &n){
	for (int i = 0; i < (int)n.GetN(); i++){
		Ptr<Node> node = n.Get(i);
		if (node->GetNodeType() == 0)
			CalculateRoute(node);
	}
}

void SetRoutingEntries(){
	// For each node.
	for (auto i = nextHop.begin(); i != nextHop.end(); i++){
		Ptr<Node> node = i->first;
		auto &table = i->second;
		for (auto j = table.begin(); j != table.end(); j++){
			// The destination node.
			Ptr<Node> dst = j->first;
			// The IP address of the dst.
			Ipv4Address dstAddr = dst->GetObject<Ipv4>()->GetAddress(1, 0).GetLocal();
			// The next hops towards the dst.
			vector<Ptr<Node> > nexts = j->second;
			for (int k = 0; k < (int)nexts.size(); k++){
				Ptr<Node> next = nexts[k];
				uint32_t interface = nbr2if[node][next].idx;
				if (node->GetNodeType() == 1)
					DynamicCast<SwitchNode>(node)->AddTableEntry(dstAddr, interface);
				else{
					node->GetObject<RdmaDriver>()->m_rdma->AddTableEntry(dstAddr, interface);
				}
			}
		}
	}
}

// take down the link between a and b, and redo the routing
void TakeDownLink(NodeContainer n, Ptr<Node> a, Ptr<Node> b){
	std::cout<<"Should never come in Takedownlink\n";
	if (!nbr2if[a][b].up)
		return;
	// take down link between a and b
	nbr2if[a][b].up = nbr2if[b][a].up = false;
	nextHop.clear();
	CalculateRoutes(n);
	// clear routing tables
	for (uint32_t i = 0; i < n.GetN(); i++){
		if (n.Get(i)->GetNodeType() == 1)
			DynamicCast<SwitchNode>(n.Get(i))->ClearTable();
		else
			n.Get(i)->GetObject<RdmaDriver>()->m_rdma->ClearTable();
	}
	DynamicCast<QbbNetDevice>(a->GetDevice(nbr2if[a][b].idx))->TakeDown();
	DynamicCast<QbbNetDevice>(b->GetDevice(nbr2if[b][a].idx))->TakeDown();
	// reset routing table
	SetRoutingEntries();

	// redistribute qp on each host
	for (uint32_t i = 0; i < n.GetN(); i++){
		if (n.Get(i)->GetNodeType() == 0)
			n.Get(i)->GetObject<RdmaDriver>()->m_rdma->RedistributeQp();
	}
}

uint64_t get_nic_rate(NodeContainer &n){
	for (uint32_t i = 0; i < n.GetN(); i++)
		if (n.Get(i)->GetNodeType() == 0)
			return DynamicCast<QbbNetDevice>(n.Get(i)->GetDevice(1))->GetDataRate().GetBitRate();
}

int main(int argc, char *argv[])
{
	clock_t begint, endt;
	begint = clock();
#ifndef PGO_TRAINING
	if (argc > 1)
#else
	if (true)
#endif
	{
		//Read the configuration file
		std::ifstream conf;
#ifndef PGO_TRAINING
		conf.open(argv[1]);
#else
		conf.open(PATH_TO_PGO_CONFIG);
#endif
		while (!conf.eof())
		{
			std::string key;
			conf >> key;

			//std::cout << conf.cur << "\n";

			if (key.compare("ENABLE_QCN") == 0)
			{
				uint32_t v;
				conf >> v;
				enable_qcn = v;
				if (enable_qcn)
					std::cout << "ENABLE_QCN\t\t\t" << "Yes" << "\n";
				else
					std::cout << "ENABLE_QCN\t\t\t" << "No" << "\n";
			}
			else if (key.compare("USE_DYNAMIC_PFC_THRESHOLD") == 0)
			{
				uint32_t v;
				conf >> v;
				use_dynamic_pfc_threshold = v;
				if (use_dynamic_pfc_threshold)
					std::cout << "USE_DYNAMIC_PFC_THRESHOLD\t" << "Yes" << "\n";
				else
					std::cout << "USE_DYNAMIC_PFC_THRESHOLD\t" << "No" << "\n";
			}
			else if (key.compare("CLAMP_TARGET_RATE") == 0)
			{
				uint32_t v;
				conf >> v;
				clamp_target_rate = v;
				if (clamp_target_rate)
					std::cout << "CLAMP_TARGET_RATE\t\t" << "Yes" << "\n";
				else
					std::cout << "CLAMP_TARGET_RATE\t\t" << "No" << "\n";
			}
			else if (key.compare("PAUSE_TIME") == 0)
			{
				double v;
				conf >> v;
				pause_time = v;
				std::cout << "PAUSE_TIME\t\t\t" << pause_time << "\n";
			}
			else if (key.compare("DATA_RATE") == 0)
			{
				std::string v;
				conf >> v;
				data_rate = v;
				std::cout << "DATA_RATE\t\t\t" << data_rate << "\n";
			}
			else if (key.compare("LINK_DELAY") == 0)
			{
				std::string v;
				conf >> v;
				link_delay = v;
				std::cout << "LINK_DELAY\t\t\t" << link_delay << "\n";
			}
			else if (key.compare("PACKET_PAYLOAD_SIZE") == 0)
			{
				uint32_t v;
				conf >> v;
				packet_payload_size = v;
				std::cout << "PACKET_PAYLOAD_SIZE\t\t" << packet_payload_size << "\n";
			}
			else if (key.compare("L2_CHUNK_SIZE") == 0)
			{
				uint32_t v;
				conf >> v;
				l2_chunk_size = v;
				std::cout << "L2_CHUNK_SIZE\t\t\t" << l2_chunk_size << "\n";
			}
			else if (key.compare("L2_ACK_INTERVAL") == 0)
			{
				uint32_t v;
				conf >> v;
				l2_ack_interval = v;
				std::cout << "L2_ACK_INTERVAL\t\t\t" << l2_ack_interval << "\n";
			}
			else if (key.compare("L2_BACK_TO_ZERO") == 0)
			{
				uint32_t v;
				conf >> v;
				l2_back_to_zero = v;
				if (l2_back_to_zero)
					std::cout << "L2_BACK_TO_ZERO\t\t\t" << "Yes" << "\n";
				else
					std::cout << "L2_BACK_TO_ZERO\t\t\t" << "No" << "\n";
			}
			else if (key.compare("TOPOLOGY_FILE") == 0)
			{
				std::string v;
				conf >> v;
				topology_file = v;
				std::cout << "TOPOLOGY_FILE\t\t\t" << topology_file << "\n";
			}
			else if (key.compare("FLOW_FILE") == 0)
			{
				std::string v;
				conf >> v;
				flow_file = v;
				std::cout << "FLOW_FILE\t\t\t" << flow_file << "\n";
			}
			else if (key.compare("TRACE_FILE") == 0)
			{
				std::string v;
				conf >> v;
				trace_file = v;
				std::cout << "TRACE_FILE\t\t\t" << trace_file << "\n";
			}
			else if (key.compare("TRACE_OUTPUT_FILE") == 0)
			{
				std::string v;
				conf >> v;
				trace_output_file = v;
				if (argc > 2)
				{
					trace_output_file = trace_output_file + std::string(argv[2]);
				}
				std::cout << "TRACE_OUTPUT_FILE\t\t" << trace_output_file << "\n";
			}
			else if (key.compare("SIMULATOR_STOP_TIME") == 0)
			{
				double v;
				conf >> v;
				simulator_stop_time = v;
				std::cout << "SIMULATOR_STOP_TIME\t\t" << simulator_stop_time << "\n";
			}
			else if (key.compare("ALPHA_RESUME_INTERVAL") == 0)
			{
				double v;
				conf >> v;
				alpha_resume_interval = v;
				std::cout << "ALPHA_RESUME_INTERVAL\t\t" << alpha_resume_interval << "\n";
			}
			else if (key.compare("RP_TIMER") == 0)
			{
				double v;
				conf >> v;
				rp_timer = v;
				std::cout << "RP_TIMER\t\t\t" << rp_timer << "\n";
			}
			else if (key.compare("EWMA_GAIN") == 0)
			{
				double v;
				conf >> v;
				ewma_gain = v;
				std::cout << "EWMA_GAIN\t\t\t" << ewma_gain << "\n";
			}
			else if (key.compare("FAST_RECOVERY_TIMES") == 0)
			{
				uint32_t v;
				conf >> v;
				fast_recovery_times = v;
				std::cout << "FAST_RECOVERY_TIMES\t\t" << fast_recovery_times << "\n";
			}
			else if (key.compare("RATE_AI") == 0)
			{
				std::string v;
				conf >> v;
				rate_ai = v;
				std::cout << "RATE_AI\t\t\t\t" << rate_ai << "\n";
			}
			else if (key.compare("RATE_HAI") == 0)
			{
				std::string v;
				conf >> v;
				rate_hai = v;
				std::cout << "RATE_HAI\t\t\t" << rate_hai << "\n";
			}
			else if (key.compare("ERROR_RATE_PER_LINK") == 0)
			{
				double v;
				conf >> v;
				error_rate_per_link = v;
				std::cout << "ERROR_RATE_PER_LINK\t\t" << error_rate_per_link << "\n";
			}
			else if (key.compare("CC_MODE") == 0){
				conf >> cc_mode;
				std::cout << "CC_MODE\t\t" << cc_mode << '\n';
			}else if (key.compare("RATE_DECREASE_INTERVAL") == 0){
				double v;
				conf >> v;
				rate_decrease_interval = v;
				std::cout << "RATE_DECREASE_INTERVAL\t\t" << rate_decrease_interval << "\n";
			}else if (key.compare("MIN_RATE") == 0){
				conf >> min_rate;
				std::cout << "MIN_RATE\t\t" << min_rate << "\n";
			}else if (key.compare("FCT_OUTPUT_FILE") == 0){
				conf >> fct_output_file;
				std::cout << "FCT_OUTPUT_FILE\t\t" << fct_output_file << '\n';
			}else if (key.compare("HAS_WIN") == 0){
				conf >> has_win;
				std::cout << "HAS_WIN\t\t" << has_win << "\n";
			}else if (key.compare("GLOBAL_T") == 0){
				conf >> global_t;
				std::cout << "GLOBAL_T\t\t" << global_t << '\n';
			}else if (key.compare("MI_THRESH") == 0){
				conf >> mi_thresh;
				std::cout << "MI_THRESH\t\t" << mi_thresh << '\n';
			}else if (key.compare("VAR_WIN") == 0){
				uint32_t v;
				conf >> v;
				var_win = v;
				std::cout << "VAR_WIN\t\t" << v << '\n';
			}else if (key.compare("FAST_REACT") == 0){
				uint32_t v;
				conf >> v;
				fast_react = v;
				std::cout << "FAST_REACT\t\t" << v << '\n';
			}else if (key.compare("U_TARGET") == 0){
				conf >> u_target;
				std::cout << "U_TARGET\t\t" << u_target << '\n';
			}else if (key.compare("INT_MULTI") == 0){
				conf >> int_multi;
				std::cout << "INT_MULTI\t\t\t\t" << int_multi << '\n';
			}else if (key.compare("RATE_BOUND") == 0){
				uint32_t v;
				conf >> v;
				rate_bound = v;
				std::cout << "RATE_BOUND\t\t" << rate_bound << '\n';
			}else if (key.compare("ACK_HIGH_PRIO") == 0){
				conf >> ack_high_prio;
				std::cout << "ACK_HIGH_PRIO\t\t" << ack_high_prio << '\n';
			}else if (key.compare("DCTCP_RATE_AI") == 0){
				conf >> dctcp_rate_ai;
				std::cout << "DCTCP_RATE_AI\t\t\t\t" << dctcp_rate_ai << "\n";
			}else if (key.compare("PFC_OUTPUT_FILE") == 0){
				conf >> pfc_output_file;
				std::cout << "PFC_OUTPUT_FILE\t\t\t\t" << pfc_output_file << '\n';
			}else if (key.compare("LINK_DOWN") == 0){
				conf >> link_down_time >> link_down_A >> link_down_B;
				std::cout << "LINK_DOWN\t\t\t\t" << link_down_time << ' '<< link_down_A << ' ' << link_down_B << '\n';
			}else if (key.compare("ENABLE_TRACE") == 0){
				conf >> enable_trace;
				std::cout << "ENABLE_TRACE\t\t\t\t" << enable_trace << '\n';
			}else if (key.compare("KMAX_MAP") == 0){
				int n_k ;
				conf >> n_k;
				std::cout << "KMAX_MAP\t\t\t\t";
				for (int i = 0; i < n_k; i++){
					uint64_t rate;
					uint32_t k;
					conf >> rate >> k;
					rate2kmax[rate] = k;
					std::cout << ' ' << rate << ' ' << k;
				}
				std::cout<<'\n';
			}else if (key.compare("KMIN_MAP") == 0){
				int n_k ;
				conf >> n_k;
				std::cout << "KMIN_MAP\t\t\t\t";
				for (int i = 0; i < n_k; i++){
					uint64_t rate;
					uint32_t k;
					conf >> rate >> k;
					rate2kmin[rate] = k;
					std::cout << ' ' << rate << ' ' << k;
				}
				std::cout<<'\n';
			}else if (key.compare("PMAX_MAP") == 0){
				int n_k ;
				conf >> n_k;
				std::cout << "PMAX_MAP\t\t\t\t";
				for (int i = 0; i < n_k; i++){
					uint64_t rate;
					double p;
					conf >> rate >> p;
					rate2pmax[rate] = p;
					std::cout << ' ' << rate << ' ' << p;
				}
				std::cout<<'\n';
			}else if (key.compare("BUFFER_SIZE") == 0){
				conf >> buffer_size;
				std::cout << "BUFFER_SIZE\t\t\t\t" << buffer_size << '\n';
			}else if (key.compare("QLEN_MON_FILE") == 0){
				conf >> qlen_mon_file;
				std::cout << "QLEN_MON_FILE\t\t\t\t" << qlen_mon_file << '\n';
			}else if (key.compare("QLEN_MON_START") == 0){
				conf >> qlen_mon_start;
				std::cout << "QLEN_MON_START\t\t\t\t" << qlen_mon_start << '\n';
			}else if (key.compare("QLEN_MON_END") == 0){
				conf >> qlen_mon_end;
				std::cout << "QLEN_MON_END\t\t\t\t" << qlen_mon_end << '\n';
			}else if (key.compare("MULTI_RATE") == 0){
				int v;
				conf >> v;
				multi_rate = v;
				std::cout << "MULTI_RATE\t\t\t\t" << multi_rate << '\n';
			}else if (key.compare("SAMPLE_FEEDBACK") == 0){
				int v;
				conf >> v;
				sample_feedback = v;
				std::cout << "SAMPLE_FEEDBACK\t\t\t\t" << sample_feedback << '\n';
			}
			else if (key.compare("WORKLOAD") == 0)
			{
				std::string v;
				conf >> v;
				if (v=="W3")
					{
						std::cout << "WORKLOAD\t\t\t" << "W3" << "\n";
						w_3 = true;
					}
				else if(v=="CONGA")
					{
						std::cout << "WORKLOAD\t\t\t" << "CONGA" << "\n";
						conga = true;
					}
					else{
						std::cout << "WORKLOAD\t\t\t" << "DCTCP" << "\n";
						w_dctcp = true;
					}
			}
			else if (key.compare("POISSON_DISTRIBUTION") == 0)
			{
				uint32_t v;
				conf >> v;
				poison_distr = v;
				if (poison_distr)
					std::cout << "POISSON_DISTRIBUTION\t\t" << "Yes" << "\n";
				else
					std::cout << "POISSON_DISTRIBUTION\t\t" << "No" << "\n";
			}
			else if (key.compare("APP_START_TIME") == 0)
			{
				double v;
				conf >> v;
				app_start_time = v;
				std::cout << "SINK_START_TIME\t\t\t" << app_start_time << "\n";
			}
			else if (key.compare("NUMBER_OF_INCAST_SOURCES") == 0)
			{
				uint32_t v;
				conf >> v;
				num_sources = v;
				std::cout << "NUMBER_OF_INCAST_SOURCES\t\t" << num_sources << "\n";
			}
			else if (key.compare("FLOW_SIZE_IN_INCAST") == 0)
			{
				uint32_t v;
				conf >> v;
				packet_size_incast = v;
				std::cout << "FLOW_SIZE_IN_INCAST\t\t" <<packet_size_incast << "\n";
			}
			else if (key.compare("CREATE_INCAST") == 0)
			{
				uint32_t v;
				conf >> v;
				create_incast = v;
				if (create_incast)
				{
					 // flow_rate -= 32 * 10 * 1000000000.0*0.05;
					 // flows_per_sec = flow_rate/(8.0*avg_flow_size);
						std::cout << "CREATE_INCAST\t\t" << "Yes" << "\n";
				}
				else
				{
					std::cout << "CREATE_INCAST\t\t" << "No" << "\n";
				}
			}
			else if (key.compare("UTILIZATION_FACTOR") == 0)
			{
				double v;
				conf >> v;
				utilization_factor = v;
				std::cout << "UTILIZATION_FACTOR\t\t\t" << utilization_factor << "\n";
			}
			else if (key.compare("APP_STOP_TIME") == 0)
			{
				double v;
				conf >> v;
				app_stop_time = v;
				std::cout << "SINK_STOP_TIME\t\t\t" << app_stop_time << "\n";
			}
			else if (key.compare("INTER_DC_TRAFFIC") == 0)
			{
				uint32_t v;
				conf >> v;
				inter_dc_traffic = v;
				if (inter_dc_traffic)
				{
					 // flow_rate -= 32 * 10 * 1000000000.0*0.05;
					 // flows_per_sec = flow_rate/(8.0*avg_flow_size);
						std::cout << "INTER_DC_TRAFFIC\t\t" << "Yes" << "\n";
				}
				else
				{
					std::cout << "INTER_DC_TRAFFIC\t\t" << "No" << "\n";
				}
			}
			else if (key.compare("QCOUNT") == 0){
				int qC ;
				conf >> qC;
				qCnt = qC;
				std::cout << "QCOUNT\t\t\t"<<qCnt<<"\n";
				std::cout<<"QCOUNT NEEDS TO BE CHANGED IN THE FOLLOWING FILES\n";
				//std::cout<<'\n';
			}
			fflush(stdout);
		}
		conf.close();
	}
	else
	{
		std::cout << "Error: require a config file\n";
		fflush(stdout);
		return 1;
	}




	bool dynamicth = use_dynamic_pfc_threshold;

	Config::SetDefault("ns3::QbbNetDevice::PauseTime", UintegerValue(pause_time));
	Config::SetDefault("ns3::QbbNetDevice::QcnEnabled", BooleanValue(enable_qcn));
	Config::SetDefault("ns3::QbbNetDevice::DynamicThreshold", BooleanValue(dynamicth));

	// set int_multi
	IntHop::multi = int_multi;
	// IntHeader::mode
	if (cc_mode == 7) // timely, use ts
		IntHeader::mode = 1;
	else if (cc_mode == 3) // hpcc, use int
		IntHeader::mode = 0;
	else // others, no extra header
		IntHeader::mode = 5;

	//SeedManager::SetSeed(time(NULL));

	// std::ifstream topof, flowf, tracef;
	// topof.open(topology_file.c_str());
	// flowf.open(flow_file.c_str());
	// tracef.open(trace_file.c_str());
	uint32_t node_num, switch_num, link_num, flow_num, trace_num;
	// topof >> node_num >> switch_num >> link_num;
	// flowf >> flow_num;
	// tracef >> trace_num;

	uint32_t tot_nodes = 154;

	node_num = 77;
	switch_num = 13;


	NodeContainer n;
	//n.Create(node_num);
	std::vector<uint32_t> node_type(tot_nodes, 0);
	for (uint32_t i = 0; i < switch_num; i++)
	{
		uint32_t sid;
		sid = i+64;
		node_type[sid] = 1;
		node_type[sid+77]=1;
	}
	for (uint32_t i = 0; i < tot_nodes; i++){
		if (node_type[i] == 0)
			n.Add(CreateObject<Node>());
		else{
			Ptr<SwitchNode> sw = CreateObject<SwitchNode>();
			n.Add(sw);
			sw->SetAttribute("EcnEnabled", BooleanValue(enable_qcn));
			//sw->SetSwitch(i);
		}
	}

	if(w_3)
{
	std::ifstream workload;
	workload.open("mix/W3.txt");

	if (!workload) {
    std::cerr << "Unable to open file W3.txt";
    exit(1);   // call system to stop
	}
	int size;
	double cdf=0.0;
	//std::cout<<"New file--------\n";
	int index = 0;
	while (workload >> size) {
		workload >> cdf;
		//std::cout<<size<<" "<<cdf<<"\n";
  		if(size>min_size)
  			{
  				enterprise_size.push_back(size);
  				enterprise_prob.push_back(cdf);
  			}
  			else{
  				enterprise_size.push_back(size);
  				enterprise_prob.push_back(0.0);
  			}
  		index++;
	}
	workload.close();
}
else if(conga)
{
	std::ifstream workload;
	workload.open("mix/w4.txt");

	if (!workload) {
    std::cerr << "Unable to open file w4.txt";
    exit(1);   // call system to stop
	}
	int size;
	double cdf=0.0;
	//std::cout<<"New file--------\n";
	int index = 0;
	while (workload >> size) {
		workload >> cdf;
		//std::cout<<size<<" "<<cdf<<"\n";
  		if(size>min_size)
  			{
  				enterprise_size.push_back(size);
  				enterprise_prob.push_back(cdf);
  			}
  			else{
  				enterprise_size.push_back(size);
  				enterprise_prob.push_back(0.0);
  			}
  		index++;
	}
	workload.close();
}
else if(w_dctcp)
{
	enterprise_size = enterprise_size_dctcp;
	enterprise_prob = enterprise_prob_dctcp;
}
for(int i=0;i<enterprise_size.size();i++)
{
	probs.insert(enterprise_prob[i]);
	size_map[enterprise_prob[i]]= enterprise_size[i];
}
double avg_flow_size=0;//in Bytes
for (int i = 1; i < enterprise_size.size(); i++) {
	avg_flow_size += (enterprise_prob[i] - enterprise_prob[i - 1]) * ((enterprise_size[i] + enterprise_size[i - 1]) / 2.0);
}	


	std::cout<<" Avg flow size "<<avg_flow_size<<"\n";

	
	std::cout<<"Utilization Factor (including incast, if there) "<<utilization_factor<<"\n";
	double flow_rate= 32 * datarate * 1000000000.0 * utilization_factor;//in Bits/s
//double flow_rate= 0.5 * 10 * 1000000000.0 * 0.6;//in Bits/s
	if(create_incast)
	{
		flow_rate -= 32*datarate*1000000000.0 *0.05;
	}
	double flows_per_sec = flow_rate / (avg_flow_size * 8.0);

	NS_LOG_INFO("Create nodes.");

	InternetStackHelper internet;
	internet.Install(n);

	//
	// Assign IP to each server
	//
	std::vector<Ipv4Address> serverAddress;
	for (uint32_t i = 0; i < tot_nodes; i++){
		if (n.Get(i)->GetNodeType() == 0){ // is server
			serverAddress.resize(i + 1);
			serverAddress[i] = Ipv4Address(0x0b000001 + ((i / 256) * 0x00010000) + ((i % 256) * 0x00000100));
		}
	}

	NS_LOG_INFO("Create channels.");

	//
	// Explicitly create the channels required by the topology.
	//

	Ptr<RateErrorModel> rem = CreateObject<RateErrorModel>();
	Ptr<UniformRandomVariable> uv = CreateObject<UniformRandomVariable>();
	rem->SetRandomVariable(uv);
	uv->SetStream(50);
	rem->SetAttribute("ErrorRate", DoubleValue(error_rate_per_link));
	rem->SetAttribute("ErrorUnit", StringValue("ERROR_UNIT_PACKET"));

	FILE *pfc_file = fopen(pfc_output_file.c_str(), "w");

	QbbHelper qbb;
	Ipv4AddressHelper ipv4;
	link_num = 208;
	for (uint32_t i = 0; i < link_num; i++)
	{
		uint32_t src, dst;
		std::string data_rate, link_delay;
		double error_rate =  0.0;
		// double error_rate;
		// topof >> src >> dst >> data_rate >> link_delay >> error_rate;
		data_rate = "100Gbps";
		link_delay = "0.001ms";
		uint32_t end_num = 128;
		if(i<104)
        {
	        if (i < 64) {
	            src=64+(i/16);
	            dst=i;
	        } else if(i<96) {
	            dst=(i-64)/8 + 64;
	            src=(i-64)%8 + 68;// Why 68? CHECK??
	        }
	        else{
	        	src = i-96+68;
	        	dst = 76;
			//data_rate = "100Gbps";
	        }
	    }
	    else{
	    	i = i -104;
	    	if (i < 64) {
	            src=64+(i/16);
	            dst=i;
	        } else if(i<96) {
	            dst=(i-64)/8 + 64;
	            src=(i-64)%8 + 68;// Why 68? CHECK??
	        }
	        else{
	        	//data_rate = "100Gbps" ;
	        	src = i-96+68;
	        	dst = 76;
	        }
	        i = i+ 104;
	        src += 77;
	        dst +=77;
	    }
		Ptr<Node> snode = n.Get(src), dnode = n.Get(dst);

		qbb.SetDeviceAttribute("DataRate", StringValue(data_rate));
		qbb.SetChannelAttribute("Delay", StringValue(link_delay));

		if (error_rate > 0)
		{
			Ptr<RateErrorModel> rem = CreateObject<RateErrorModel>();
			Ptr<UniformRandomVariable> uv = CreateObject<UniformRandomVariable>();
			rem->SetRandomVariable(uv);
			uv->SetStream(50);
			rem->SetAttribute("ErrorRate", DoubleValue(error_rate));
			rem->SetAttribute("ErrorUnit", StringValue("ERROR_UNIT_PACKET"));
			qbb.SetDeviceAttribute("ReceiveErrorModel", PointerValue(rem));
		}
		else
		{
			qbb.SetDeviceAttribute("ReceiveErrorModel", PointerValue(rem));
		}

		fflush(stdout);

		// Assigne server IP
		// Note: this should be before the automatic assignment below (ipv4.Assign(d)),
		// because we want our IP to be the primary IP (first in the IP address list),
		// so that the global routing is based on our IP
		NetDeviceContainer d = qbb.Install(snode, dnode);
		if (snode->GetNodeType() == 0){
			Ptr<Ipv4> ipv4 = snode->GetObject<Ipv4>();
			ipv4->AddInterface(d.Get(0));
			ipv4->AddAddress(1, Ipv4InterfaceAddress(serverAddress[src], Ipv4Mask(0xff000000)));
		}
		if (dnode->GetNodeType() == 0){
			Ptr<Ipv4> ipv4 = dnode->GetObject<Ipv4>();
			ipv4->AddInterface(d.Get(1));
			ipv4->AddAddress(1, Ipv4InterfaceAddress(serverAddress[dst], Ipv4Mask(0xff000000)));
		}

		// used to create a graph of the topology
		nbr2if[snode][dnode].idx = DynamicCast<QbbNetDevice>(d.Get(0))->GetIfIndex();
		nbr2if[snode][dnode].up = true;
		nbr2if[snode][dnode].delay = DynamicCast<QbbChannel>(DynamicCast<QbbNetDevice>(d.Get(0))->GetChannel())->GetDelay().GetTimeStep();
		nbr2if[snode][dnode].bw = DynamicCast<QbbNetDevice>(d.Get(0))->GetDataRate().GetBitRate();
		nbr2if[dnode][snode].idx = DynamicCast<QbbNetDevice>(d.Get(1))->GetIfIndex();
		nbr2if[dnode][snode].up = true;
		nbr2if[dnode][snode].delay = DynamicCast<QbbChannel>(DynamicCast<QbbNetDevice>(d.Get(1))->GetChannel())->GetDelay().GetTimeStep();
		nbr2if[dnode][snode].bw = DynamicCast<QbbNetDevice>(d.Get(1))->GetDataRate().GetBitRate();

		// This is just to set up the connectivity between nodes. The IP addresses are useless
		char ipstring[16];
		sprintf(ipstring, "10.%d.%d.0", i / 254 + 1, i % 254 + 1);
		ipv4.SetBase(ipstring, "255.255.255.0");
		ipv4.Assign(d);

		// setup PFC trace
		DynamicCast<QbbNetDevice>(d.Get(0))->TraceConnectWithoutContext("QbbPfc", MakeBoundCallback (&get_pfc, pfc_file, DynamicCast<QbbNetDevice>(d.Get(0))));
		DynamicCast<QbbNetDevice>(d.Get(1))->TraceConnectWithoutContext("QbbPfc", MakeBoundCallback (&get_pfc, pfc_file, DynamicCast<QbbNetDevice>(d.Get(1))));
	}

	//INTER DC LINK
	uint32_t src, dst, i;
	src = 76;
	dst = 153;
	i = 209;

	Ptr<Node> snode = n.Get(src), dnode = n.Get(dst);

	qbb.SetDeviceAttribute("DataRate", StringValue("200Gbps"));
	qbb.SetChannelAttribute("Delay", StringValue("0.2ms"));

	qbb.SetDeviceAttribute("ReceiveErrorModel", PointerValue(rem));

	fflush(stdout);

		// Assigne server IP
		// Note: this should be before the automatic assignment below (ipv4.Assign(d)),
		// because we want our IP to be the primary IP (first in the IP address list),
		// so that the global routing is based on our IP
		NetDeviceContainer d = qbb.Install(snode, dnode);
		if (snode->GetNodeType() == 0){
			Ptr<Ipv4> ipv4 = snode->GetObject<Ipv4>();
			ipv4->AddInterface(d.Get(0));
			ipv4->AddAddress(1, Ipv4InterfaceAddress(serverAddress[src], Ipv4Mask(0xff000000)));
		}
		if (dnode->GetNodeType() == 0){
			Ptr<Ipv4> ipv4 = dnode->GetObject<Ipv4>();
			ipv4->AddInterface(d.Get(1));
			ipv4->AddAddress(1, Ipv4InterfaceAddress(serverAddress[dst], Ipv4Mask(0xff000000)));
		}

		// used to create a graph of the topology
		nbr2if[snode][dnode].idx = DynamicCast<QbbNetDevice>(d.Get(0))->GetIfIndex();
		nbr2if[snode][dnode].up = true;
		nbr2if[snode][dnode].delay = DynamicCast<QbbChannel>(DynamicCast<QbbNetDevice>(d.Get(0))->GetChannel())->GetDelay().GetTimeStep();
		nbr2if[snode][dnode].bw = DynamicCast<QbbNetDevice>(d.Get(0))->GetDataRate().GetBitRate();
		nbr2if[dnode][snode].idx = DynamicCast<QbbNetDevice>(d.Get(1))->GetIfIndex();
		nbr2if[dnode][snode].up = true;
		nbr2if[dnode][snode].delay = DynamicCast<QbbChannel>(DynamicCast<QbbNetDevice>(d.Get(1))->GetChannel())->GetDelay().GetTimeStep();
		nbr2if[dnode][snode].bw = DynamicCast<QbbNetDevice>(d.Get(1))->GetDataRate().GetBitRate();

		// This is just to set up the connectivity between nodes. The IP addresses are useless
		char ipstring[16];
		sprintf(ipstring, "10.%d.%d.0", i / 254 + 1, i % 254 + 1);
		ipv4.SetBase(ipstring, "255.255.255.0");
		ipv4.Assign(d);

		// setup PFC trace
		DynamicCast<QbbNetDevice>(d.Get(0))->TraceConnectWithoutContext("QbbPfc", MakeBoundCallback (&get_pfc, pfc_file, DynamicCast<QbbNetDevice>(d.Get(0))));
		DynamicCast<QbbNetDevice>(d.Get(1))->TraceConnectWithoutContext("QbbPfc", MakeBoundCallback (&get_pfc, pfc_file, DynamicCast<QbbNetDevice>(d.Get(1))));



	/////////////////////
	nic_rate = get_nic_rate(n);

	// config switch
	for (uint32_t i = 0; i < tot_nodes; i++){
		if (n.Get(i)->GetNodeType() == 1){ // is switch
			Ptr<SwitchNode> sw = DynamicCast<SwitchNode>(n.Get(i));
			uint32_t shift = 3; // by default 1/8
			for (uint32_t j = 1; j < sw->GetNDevices(); j++){
				Ptr<QbbNetDevice> dev = DynamicCast<QbbNetDevice>(sw->GetDevice(j));
				// set ecn
				uint64_t rate = dev->GetDataRate().GetBitRate();
				NS_ASSERT_MSG(rate2kmin.find(rate) != rate2kmin.end(), "must set kmin for each link speed");
				NS_ASSERT_MSG(rate2kmax.find(rate) != rate2kmax.end(), "must set kmax for each link speed");
				NS_ASSERT_MSG(rate2pmax.find(rate) != rate2pmax.end(), "must set pmax for each link speed");
				sw->m_mmu->ConfigEcn(j, rate2kmin[rate], rate2kmax[rate], rate2pmax[rate]);
				// set pfc
				uint64_t delay = DynamicCast<QbbChannel>(dev->GetChannel())->GetDelay().GetTimeStep();
				uint32_t headroom = rate * delay / 8 / 1000000000 * 3;
				sw->m_mmu->ConfigHdrm(j, headroom);

				// set pfc alpha, proportional to link bw
				sw->m_mmu->pfc_a_shift[j] = shift;
				while (rate > nic_rate && sw->m_mmu->pfc_a_shift[j] > 0){
					sw->m_mmu->pfc_a_shift[j]--;
					rate /= 2;
				}
			}
			sw->m_mmu->ConfigNPort(sw->GetNDevices()-1);
			sw->m_mmu->ConfigBufferSize(buffer_size* 1024 * 1024);
			sw->m_mmu->node_id = sw->GetId();
			sw->m_mmu->node_id = sw->GetId();
			sw->m_mmu->SetSwitch(i);
			if((i+1)%77==0)
			{
				sw->m_mmu->SetGateway();
			}
		}

	}

	#if ENABLE_QP
	FILE *fct_output = fopen(fct_output_file.c_str(), "w");
	//
	// install RDMA driver
	//
	for (uint32_t i = 0; i < tot_nodes; i++){
		if (n.Get(i)->GetNodeType() == 0){ // is server
			// create RdmaHw
			Ptr<RdmaHw> rdmaHw = CreateObject<RdmaHw>();
			rdmaHw->SetAttribute("ClampTargetRate", BooleanValue(clamp_target_rate));
			rdmaHw->SetAttribute("AlphaResumInterval", DoubleValue(alpha_resume_interval));
			rdmaHw->SetAttribute("RPTimer", DoubleValue(rp_timer));
			rdmaHw->SetAttribute("FastRecoveryTimes", UintegerValue(fast_recovery_times));
			rdmaHw->SetAttribute("EwmaGain", DoubleValue(ewma_gain));
			rdmaHw->SetAttribute("RateAI", DataRateValue(DataRate(rate_ai)));
			rdmaHw->SetAttribute("RateHAI", DataRateValue(DataRate(rate_hai)));
			rdmaHw->SetAttribute("L2BackToZero", BooleanValue(l2_back_to_zero));
			rdmaHw->SetAttribute("L2ChunkSize", UintegerValue(l2_chunk_size));
			rdmaHw->SetAttribute("L2AckInterval", UintegerValue(l2_ack_interval));
			rdmaHw->SetAttribute("CcMode", UintegerValue(cc_mode));
			rdmaHw->SetAttribute("RateDecreaseInterval", DoubleValue(rate_decrease_interval));
			rdmaHw->SetAttribute("MinRate", DataRateValue(DataRate(min_rate)));
			rdmaHw->SetAttribute("Mtu", UintegerValue(packet_payload_size));
			rdmaHw->SetAttribute("MiThresh", UintegerValue(mi_thresh));
			rdmaHw->SetAttribute("VarWin", BooleanValue(var_win));
			rdmaHw->SetAttribute("FastReact", BooleanValue(fast_react));
			rdmaHw->SetAttribute("MultiRate", BooleanValue(multi_rate));
			rdmaHw->SetAttribute("SampleFeedback", BooleanValue(sample_feedback));
			rdmaHw->SetAttribute("TargetUtil", DoubleValue(u_target));
			rdmaHw->SetAttribute("RateBound", BooleanValue(rate_bound));
			rdmaHw->SetAttribute("DctcpRateAI", DataRateValue(DataRate(dctcp_rate_ai)));
			// create and install RdmaDriver
			Ptr<RdmaDriver> rdma = CreateObject<RdmaDriver>();
			Ptr<Node> node = n.Get(i);
			rdma->SetNode(node);
			rdma->SetRdmaHw(rdmaHw);

			node->AggregateObject (rdma);
			rdma->Init();
			rdma->TraceConnectWithoutContext("QpComplete", MakeBoundCallback (qp_finish, fct_output));
		}
	}
	#endif

	// set ACK priority on hosts
	if (ack_high_prio)
		RdmaEgressQueue::ack_q_idx = 0;
	else
		RdmaEgressQueue::ack_q_idx = 3;

	//
	// setup switch CC
	//
	for (uint32_t i = 0; i < tot_nodes; i++){
		if (n.Get(i)->GetNodeType() == 1){ // switch
			Ptr<SwitchNode> sw = DynamicCast<SwitchNode>(n.Get(i));
			sw->SetAttribute("CcMode", UintegerValue(cc_mode));
		}
	}

	// setup routing
	CalculateRoutes(n);
	SetRoutingEntries();

	//
	// get BDP and delay
	//
	maxRtt = maxBdp = 0;
	for (uint32_t i = 0; i < tot_nodes; i++){
		if (n.Get(i)->GetNodeType() != 0)
			continue;
		for (uint32_t j = i+1; j < tot_nodes; j++){
			if (n.Get(j)->GetNodeType() != 0)
				continue;
			uint64_t delay = pairDelay[n.Get(i)][n.Get(j)];
			uint64_t txDelay = pairTxDelay[n.Get(i)][n.Get(j)];
			uint64_t rtt = delay * 2 + txDelay;
			uint64_t bw = pairBw[n.Get(i)][n.Get(j)];
			uint64_t bdp = rtt * bw / 1000000000/8; 
			pairBdp[n.Get(i)][n.Get(j)] = bdp;
			pairRtt[n.Get(i)][n.Get(j)] = rtt;
			if (bdp > maxBdp)
				maxBdp = bdp;
			if (rtt > maxRtt)
				maxRtt = rtt;
		}
	}
	printf("maxRtt=%lu maxBdp=%lu\n", maxRtt, maxBdp);

	//
	// add trace
	//

	// NodeContainer trace_nodes;
	// for (uint32_t i = 0; i < trace_num; i++)
	// {
	// 	uint32_t nid;
	// 	tracef >> nid;
	// 	if (nid >= n.GetN()){
	// 		continue;
	// 	}
	// 	trace_nodes = NodeContainer(trace_nodes, n.Get(nid));
	// }

	// FILE *trace_output = fopen(trace_output_file.c_str(), "w");
	// if (enable_trace)
	// 	qbb.EnableTracing(trace_output, trace_nodes);

	// // dump link speed to trace file
	// {
	// 	SimSetting sim_setting;
	// 	for (auto i: nbr2if){
	// 		for (auto j : i.second){
	// 			uint16_t node = i.first->GetId();
	// 			uint8_t intf = j.second.idx;
	// 			uint64_t bps = DynamicCast<QbbNetDevice>(i.first->GetDevice(j.second.idx))->GetDataRate().GetBitRate();
	// 			sim_setting.port_speed[node][intf] = bps;
	// 		}
	// 	}
	// 	sim_setting.win = maxBdp;
	// 	sim_setting.Serialize(trace_output);
	// }
	uint32_t packetSize = packet_payload_size;
	//Time interPacketInterval = Seconds(0.0000005 / 2);
	double interarrival = (packetSize*8)/(1000000000.0*datarate);
	//std::cout<<interarrival<<"in\n";
	Time interPacketInterval = Seconds(interarrival); //800ns assuming 10Gbps link and 1000 Byte Payload
    std::mt19937 gen(5489U); //Same Seed for all the Simlulations
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::exponential_distribution<double> exp_dis(flows_per_sec);
    double m = -(sigma*sigma/2)-log(flows_per_sec);
    std::lognormal_distribution<double> lognormal_dis(m,sigma);
	double start_time=app_start_time;
    //for (uint32_t i = 0; i < flow_num; i++)
    uint32_t flownum=0;
    
    uint32_t num_of_incasts = (app_stop_time-app_start_time)*32*datarate*1000000000.0*0.05/(packet_size_incast*packetSize*num_sources*8.0);
    uint32_t incasts_done = 0;
    double incast_interval = (app_stop_time-app_start_time)/num_of_incasts;
    if(create_incast)
    {
    	std::cout<<"Incast info- "<<"Packet size "<<packet_size_incast<<"Num of sources "<<num_sources<<"\n";
    }

	Ipv4GlobalRoutingHelper::PopulateRoutingTables();

	NS_LOG_INFO("Create Applications.");

	//Time interPacketInterval = Seconds(0.0000005 / 2);

	// maintain port number for each host
	std::unordered_map<uint32_t, uint16_t> portNumder;
	std::unordered_map<uint32_t, uint16_t> dportNumder;
	for (uint32_t i = 0; i < tot_nodes; i++){
		if (n.Get(i)->GetNodeType() == 0)
			portNumder[i] = 10000; // each host use port number from 10000
	}
	// for (uint32_t i = 0; i < node_num; i++){
	// 	if (n.Get(i)->GetNodeType() == 0)
	// 		portNumder[i] = 10000; // each host use port number from 10000
	// }

	// for (uint32_t i = 0; i < flow_num; i++)
	// {
	// 	uint32_t src, dst, pg, maxPacketCount, port, dport;
	// 	double start_time, stop_time;
	// 	flowf >> src >> dst >> pg >> dport >> maxPacketCount >> start_time;
	// 	NS_ASSERT(n.Get(src)->GetNodeType() == 0 && n.Get(dst)->GetNodeType() == 0);
	// 	port = portNumder[src]++; // get a new port number 
	// 	RdmaClientHelper clientHelper(pg, serverAddress[src], serverAddress[dst], port, dport, maxPacketCount, has_win?(global_t==1?maxBdp:pairBdp[n.Get(src)][n.Get(dst)]):0, global_t==1?maxRtt:pairRtt[n.Get(src)][n.Get(dst)]);
	// 	ApplicationContainer appCon = clientHelper.Install(n.Get(src));
	// 	appCon.Start(Seconds(start_time));
	// 	appCon.Stop(Seconds(stop_time));
	// }

	while (start_time<app_stop_time)
    {
        flownum += 1;
		//uint32_t src, dst, pg, maxPacketCount, port;
		uint32_t src, dst, pg, maxPacketCount, port, dport;
        src = uint32_t(dis(gen) * 64);
        //src = 1;
        while (true) {
            dst = uint32_t(dis(gen) * 64);
            if (dst != src) break;
        }
        NS_ASSERT(dst < 144);
        NS_ASSERT(src < 144);
       	pg = 3+(dis(gen)*(qCnt - 8));
       	//} //priority between 1 & 125, 0, 126 & 127 reserved
        double flow_size_helper = dis(gen);
        int start_ind=0, end_ind=enterprise_size.size();
        auto higher_prob = probs.upper_bound(flow_size_helper);
        double higher_probability = *(higher_prob);
        higher_prob--;
        double lower_probability = *(higher_prob);
        int higher_size = size_map[higher_probability];
        int size = size_map[lower_probability];
        //std::cout<<"Hs "<<higher_size<<"size "<<size<<"\n";
    	int flow_size = size + (higher_size-size)*((higher_probability - flow_size_helper)/(higher_probability - lower_probability));
    	//std::cout<<flow_size<<"\n";
    	int flow_packet_size;
    	if(flow_size <= packetSize)
    	{
    		flow_packet_size = flow_size;
    		maxPacketCount = 1;
    	}
    	else{
    		flow_packet_size = packetSize;
    		maxPacketCount = int(flow_size/packetSize)+1;
    	}
	
   
        if(flownum%1000==0) std::cout<<"New Flow Created Src "<<src<<" Dst "<<dst<<" FlowNum "<<flownum<<"Packets "<<maxPacketCount<<" Priority "<<pg<<"\n";
		NS_ASSERT(n.Get(src)->GetNodeType() == 0 && n.Get(dst)->GetNodeType() == 0);
		port = portNumder[src]++; // get a new port number 
		port = port%65000;
		if(port<10000) port+= 30000;
		RdmaClientHelper clientHelper(pg, serverAddress[src], serverAddress[dst], port, 40000+dst, maxPacketCount*flow_packet_size, has_win?(global_t==1?maxBdp:pairBdp[n.Get(src)][n.Get(dst)]):0, global_t==1?maxRtt:pairRtt[n.Get(src)][n.Get(dst)]);
		ApplicationContainer appCon = clientHelper.Install(n.Get(src));
		appCon.Start(Seconds(start_time));
		appCon.Stop(Seconds(app_stop_time));
		if(create_incast)
        {
        	if((uint32_t)((start_time- app_start_time)/(incast_interval))>incasts_done)
        	{
        		uint32_t dst2 = dst;
        		while (true)
	        	{
	            		dst = uint32_t(dis(gen) * 64) ;
	           		 	if (dst != dst2) break;
	       		}
	        	incasts_done++;
	        	//std::cout<<"Incast done "<<incasts_done<<" \n";
	        	//dst = uint32_t(dis(gen) * 64);
	        	for(int k = 0; k < num_sources;k++)
	        	{
	        		while (true)
	        		 {
	            		src = uint32_t(dis(gen) * 64) ;
	           		 	if (dst != src) break;
	       			 }
	       			 maxPacketCount= packet_size_incast;
	       			pg = 3+(dis(gen)*(qCnt - 8));
	       			
	       			NS_ASSERT(n.Get(src)->GetNodeType() == 0 && n.Get(dst)->GetNodeType() == 0);
		port = portNumder[src]++; // get a new port number 
		port = port%65000;
		if(port<10000) port+= 30000;
		RdmaClientHelper clientHelper(pg, serverAddress[src], serverAddress[dst], port, 40000+dst, maxPacketCount*packetSize, has_win?(global_t==1?maxBdp:pairBdp[n.Get(src)][n.Get(dst)]):0, global_t==1?maxRtt:pairRtt[n.Get(src)][n.Get(dst)]);
		ApplicationContainer appCon = clientHelper.Install(n.Get(src));
		appCon.Start(Seconds(start_time));
		appCon.Stop(Seconds(app_stop_time));
	        	}
	        }
        }
        double x;
        if(poison_distr)
        {
        	x = exp_dis(gen);
        	start_time += x;
        }
        else{
        	x = lognormal_dis(gen) ;
        	start_time += x;
        }
        //std::cout<<x<<"\n";
        
	}

	incasts_done = 0;
	start_time = app_start_time;

	while (start_time<app_stop_time)
    {
        flownum += 1;
		//uint32_t src, dst, pg, maxPacketCount, port;
		uint32_t src, dst, pg, maxPacketCount, port, dport;
        src = uint32_t(dis(gen) * 64)+77;
        //src = 1;
        while (true) {
            dst = uint32_t(dis(gen) * 64)+77;
            if (dst != src) break;
        }
        // NS_ASSERT(dst < 144);
        // NS_ASSERT(src < 144);
       	pg = 3+(dis(gen)*(qCnt - 8));
       	//} //priority between 1 & 125, 0, 126 & 127 reserved
        double flow_size_helper = dis(gen);
        int start_ind=0, end_ind=enterprise_size.size();
        auto higher_prob = probs.upper_bound(flow_size_helper);
        double higher_probability = *(higher_prob);
        higher_prob--;
        double lower_probability = *(higher_prob);
        int higher_size = size_map[higher_probability];
        int size = size_map[lower_probability];
        //std::cout<<"Hs "<<higher_size<<"size "<<size<<"\n";
    	int flow_size = size + (higher_size-size)*((higher_probability - flow_size_helper)/(higher_probability - lower_probability));
    	//std::cout<<flow_size<<"\n";
    	int flow_packet_size;
    	if(flow_size <= packetSize)
    	{
    		flow_packet_size = flow_size;
    		maxPacketCount = 1;
    	}
    	else{
    		flow_packet_size = packetSize;
    		maxPacketCount = int(flow_size/packetSize)+1;
    	}
	
   
        if(flownum%1000==0) std::cout<<"New Flow Created Src "<<src<<" Dst "<<dst<<" FlowNum "<<flownum<<"Packets "<<maxPacketCount<<" Priority "<<pg<<"\n";
		NS_ASSERT(n.Get(src)->GetNodeType() == 0 && n.Get(dst)->GetNodeType() == 0);
		port = portNumder[src]++; // get a new port number 
		port = port%65000;
		if(port<10000) port+= 30000;
		RdmaClientHelper clientHelper(pg, serverAddress[src], serverAddress[dst], port, 40000+dst, maxPacketCount*flow_packet_size, has_win?(global_t==1?maxBdp:pairBdp[n.Get(src)][n.Get(dst)]):0, global_t==1?maxRtt:pairRtt[n.Get(src)][n.Get(dst)]);
		ApplicationContainer appCon = clientHelper.Install(n.Get(src));
		appCon.Start(Seconds(start_time));
		appCon.Stop(Seconds(app_stop_time));
		if(create_incast)
        {
        	if((uint32_t)((start_time- app_start_time)/(incast_interval))>incasts_done)
        	{
        		uint32_t dst2 = dst;
        		while (true)
	        	{
	            		dst = uint32_t(dis(gen) * 64) +77;
	           		 	if (dst != dst2) break;
	       		}
	        	incasts_done++;
	        	//std::cout<<"Incast done "<<incasts_done<<" \n";
	        	//dst = uint32_t(dis(gen) * 64);
	        	for(int k = 0; k < num_sources;k++)
	        	{
	        		while (true)
	        		 {
	            		src = uint32_t(dis(gen) * 64)+77 ;
	           		 	if (dst != src) break;
	       			 }
	       			 maxPacketCount= packet_size_incast;
	       			pg = 3+(dis(gen)*(qCnt - 8));
	       			
	       			NS_ASSERT(n.Get(src)->GetNodeType() == 0 && n.Get(dst)->GetNodeType() == 0);
		port = portNumder[src]++; // get a new port number 
		port = port%65000;
		if(port<10000) port+= 30000;
		RdmaClientHelper clientHelper(pg, serverAddress[src], serverAddress[dst], port, 40000+dst, maxPacketCount*packetSize, has_win?(global_t==1?maxBdp:pairBdp[n.Get(src)][n.Get(dst)]):0, global_t==1?maxRtt:pairRtt[n.Get(src)][n.Get(dst)]);
		ApplicationContainer appCon = clientHelper.Install(n.Get(src));
		appCon.Start(Seconds(start_time));
		appCon.Stop(Seconds(app_stop_time));
	        	}
	        }
        }
        double x;
        if(poison_distr)
        {
        	x = exp_dis(gen);
        	start_time += x;
        }
        else{
        	x = lognormal_dis(gen) ;
        	start_time += x;
        }
        //std::cout<<x<<"\n";
        
	}

	//FLOWS IN SECOND DC


	/////INTER DC FLOWS

	start_time = app_start_time;
	if(inter_dc_traffic)
	{
		start_time = app_start_time;
		int tot_inter = 10;
		std::set<int> nodes;
	    //std::cout<<"StartTime "<<start_time<<" App Start "<<app_start_time<<"\n";
	    for(int inter_flows = 0;inter_flows<tot_inter;inter_flows++)
	    {
	        flownum += 1;
			//uint32_t src, dst, pg, maxPacketCount, port;
			uint32_t src, dst, pg, maxPacketCount, port, dport;
	        while (true) {
	            src = uint32_t(dis(gen) * 64);
	            if (nodes.find(src)==nodes.end()) {
	            	nodes.insert(src);
	            	break;
	            }
	        }
	        //src = 1;
	        while (true) {
	            dst = uint32_t(dis(gen) * 64);
	            if (dst != src && nodes.find(dst)==nodes.end()) {
	            	nodes.insert(dst);
	            	break;
	            }
	        }
	        // NS_ASSERT(dst < 64);
	        // NS_ASSERT(src < 64);
	        if(inter_flows< tot_inter/2)
	        {
	        	src += 77;
	        }
	        else{
	        	dst += 77;
	        }
	        // NS_ASSERT(dst < 144);
	        // NS_ASSERT(src < 144);
	       	pg = 3+(dis(gen)*(qCnt - 8));
	       	//} //priority between 1 & 125, 0, 126 & 127 reserved
	        int flow_packet_size = packetSize;
	    	maxPacketCount = 10000000;
		
	   
	        if(flownum%1000==0) std::cout<<"New Flow Created Src "<<src<<" Dst "<<dst<<" FlowNum "<<flownum<<"Packets "<<maxPacketCount<<" Priority "<<pg<<"\n";
			NS_ASSERT(n.Get(src)->GetNodeType() == 0 && n.Get(dst)->GetNodeType() == 0);
			port = portNumder[src]++; // get a new port number 
			port = port%65000;
			if(port<10000) port+= 30000;
			RdmaClientHelper clientHelper(pg, serverAddress[src], serverAddress[dst], port, 40000+dst, maxPacketCount*flow_packet_size, has_win?(global_t==1?maxBdp:pairBdp[n.Get(src)][n.Get(dst)]):0, global_t==1?maxRtt:pairRtt[n.Get(src)][n.Get(dst)]);
			ApplicationContainer appCon = clientHelper.Install(n.Get(src));
			appCon.Start(Seconds(start_time));
			appCon.Stop(Seconds(app_stop_time));
	        //std::cout<<x<<"\n";
	        
		}
	}

	////////////

	std::cout<<"Number of flows "<<flownum<<"\n";
	if(create_incast)
	{
		std::cout<<"Number of incasts done "<<incasts_done<<"\n";
	}


	// topof.close();
	// flowf.close();
	// tracef.close();

	// schedule link down
	if (link_down_time > 0){
		Simulator::Schedule(Seconds(2) + MicroSeconds(link_down_time), &TakeDownLink, n, n.Get(link_down_A), n.Get(link_down_B));
	}

	// schedule buffer monitor
	FILE* qlen_output = fopen(qlen_mon_file.c_str(), "w");
	Simulator::Schedule(NanoSeconds(qlen_mon_start), &monitor_buffer, qlen_output, &n);

	//
	// Now, do the actual simulation.
	//
	std::cout << "Running Simulation.\n";
	fflush(stdout);
	NS_LOG_INFO("Run Simulation.");
	Simulator::Stop(Seconds(simulator_stop_time));
	Simulator::Run();
	Simulator::Destroy();
	NS_LOG_INFO("Done.");
	//fclose(trace_output);

	endt = clock();
	std::cout << (double)(endt - begint) / CLOCKS_PER_SEC << "\n";

}
