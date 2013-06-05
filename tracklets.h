#include <iostream>
#include <ostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <stdlib.h> 
#include <cmath>
#include <set>
#include <algorithm> 
using namespace std;

struct trackLet {
        double m_dx;
        double m_dy;

        double m_ex;
        double m_ey;

       int trackLetNum ;
       vector<int> hits;
};
struct trackCand{
	int trackNum;
	vector<int> trackLets;
};
struct track {
       // double m_x0;
        //double m_dx;
       // double m_y0;
       // double m_dy;

       // double m_ex;
       // double m_ey;

       int trackHitsNum;
       set<int> Hits;    
};
 int num_events;
 int max_hits;
 int sens_num;

typedef struct {
        int startPosition;
        int hitsNum;
        double z;
} sensorInfo;

void readFile(string filename, char*& input, bool do_malloc);
void searchTrackLets(int event_no, vector<trackLet>& vector_trackLets);
void setTrackLets(trackLet *tr, int hit0offset, int hit1offset);
void setTrack(track *trk, int hit0, int hit1);
void trackLetsSelection(vector<trackLet>& vector_trackLets, vector<track>& trackCandidates);
void findtracklet(vector<trackLet>::iterator& first,vector<trackLet>::iterator& last, int val);
// findClusterdxdy
// findClusterexey
// vector<vector<trackLet> > trackCandidates_vector;
// findTracks
// vector<vector<trackLet> > tracks_vector;
vector<vector<trackLet> > track_lets_vector;
vector<vector<track> > trackCand_vector;
template<class T>
void printResultTracks(vector<T> tracks, int event_no, string track_folder);

template<class T>
void clean_vector(vector<vector<T> >& _tracks_vector);
template<class T>
void clean_vector2(vector<T>& _tracks_vector);

int* numHitsEvents;
int* sensor_hitStarts;
int* sensor_hitsNums;
double* sensor_Zs;
int* hit_ids;
int* track_ids;
int* hit_sensorNums;
int* hit_isUseds;
double* hit_Xs;
double* hit_Ys;
double* hit_Zs;
double* hit_Ws;
template <class T>
static std::string toString(T t){
        std::stringstream ss;
        std::string s;
        ss << t;
        ss >> s;
        return s;
}
/*
int main(int argc, char* argv[]){

	char c1, c2, *optarg;
        int num_events = 200;
        int sens_num = 48;
        int max_hits = 3792;

        int num_events_to_process = 4;


	string c = "input.dump";
        char* input;

       // Read file
        readFile(c, input, true);

	int dataOffset = max_hits * num_events;

        numHitsEvents = (int*) &input[0];
        sensor_hitStarts = (int*) (numHitsEvents + num_events);
        sensor_hitsNums = (int*) (sensor_hitStarts + num_events * sens_num);
        sensor_Zs = (double*) (sensor_hitsNums + num_events * sens_num);

        hit_ids = (int*) (sensor_Zs + num_events * sens_num);
        track_ids = hit_ids +  dataOffset;
        hit_sensorNums = track_ids + dataOffset;
        hit_isUseds = hit_sensorNums + dataOffset;
        hit_Xs = (double*) (hit_isUseds + dataOffset);
        hit_Ys = hit_Xs + dataOffset;
        hit_Zs = hit_Ys + dataOffset;
        hit_Ws = hit_Zs + dataOffset;
                

	int* numHitsEvents;
	numHitsEvents = (int*) &input[0];	
	int numerohit = *numHitsEvents;

        //float Whatis[200];    
 
 track_lets_vector.clear();
 trackCand_vector.clear();
        for (int i=0; i<num_events_to_process; ++i){
		vector<trackLet> t1;
		vector<track> t2;
                track_lets_vector.push_back(t1);
		trackCand_vector.push_back(t2);
                
        }

readFile(c, input, false);
clean_vector(track_lets_vector);
clean_vector(trackCand_vector);


	for (int i=1; i<=num_events_to_process; i++){

        
      	 searchTrackLets(i-1, track_lets_vector[i-1]);
	 trackLetsSelection(track_lets_vector[i-1], trackCand_vector[i-1]);

cout << "Writing results..." << endl;
//for (int i=0; i<num_events_to_process; i++)
        printResultTracks(trackCand_vector[i-1], i-1, "tracks");


}

 //               }

  }//end




void readFile(string filename, char*& input, bool do_malloc ){
        int size;

       // Give me them datas!!11!
         ifstream infile (filename.c_str(), ifstream::binary);
       // get size of file
         infile.seekg(0, ifstream::end);
         size = infile.tellg();
         infile.seekg(0);
        // read content of infile with pointers
     if (do_malloc)
         input = (char*) malloc(size);
         infile.read (input, size);
         infile.close();
  }

void searchTrackLets(int eventRealId, vector<trackLet>& trackLets_vector) {
//my code
	 int sens_num = 48 ;
	 int max_hits = 3792;
	 int num_events = 200;

	int firstSensor = 0;
	int lastSensor  = sens_num-3;
        
        int eventId = eventRealId % num_events;
     

        trackLet m_trackLet;
	m_trackLet.trackLetNum =0;

        int event_sensor_displ = eventId * sens_num;
        int event_hit_displ = eventId * max_hits;

        int sens0, sens1, first1, hit0_no, hit1_no;
        sensorInfo sensor0, sensor1, extra_sensor;
        double dxMax, dyMax;

	for ( sens0 = firstSensor; lastSensor >= sens0; sens0 += 1 ) {

	 sens1 = sens0 + 2;

	 sensor0.startPosition = sensor_hitStarts[event_sensor_displ + sens0];
         sensor0.hitsNum        = sensor_hitsNums[event_sensor_displ + sens0];
         sensor0.z = sensor_Zs[event_sensor_displ + sens0];

         sensor1.startPosition = sensor_hitStarts[event_sensor_displ + sens1];
         sensor1.hitsNum        = sensor_hitsNums[event_sensor_displ + sens1];
         sensor1.z = sensor_Zs[event_sensor_displ + sens1];

	 dxMax = 0.4 * fabs( sensor1.z - sensor0.z );
         dyMax = 0.3 * fabs( sensor1.z - sensor0.z );
         first1 = 0;

         for (hit0_no = 0; hit0_no < sensor0.hitsNum; hit0_no++){

	  int hit0_offset = event_hit_displ + sensor0.startPosition + hit0_no;
	  double x0  = hit_Xs[hit0_offset];
          double y0  = hit_Ys[hit0_offset];
          double z0  = hit_Zs[hit0_offset];

	for (hit1_no = 0; hit1_no < sensor1.hitsNum; hit1_no++) {

                        int hit1_offset = event_hit_displ + sensor1.startPosition + hit1_no;

                        double x1  = hit_Xs[hit1_offset];
			double y1  = hit_Ys[hit1_offset];
	                double z1  = hit_Zs[hit1_offset];

	         	double dx =  (x1 -x0)/(z1 - z0) ;
			double dy =  (y1 -y0)/(z1 - z0) ;
	if( dx > 0.4 || dx < -0.4){
                 continue;
                }
        if( dy > 0.3 ||  dy < -0.3){
                 continue;
                }

	// m_trackLet.hits.clear();
	
	setTrackLets(&m_trackLet, hit0_offset, hit1_offset);

//	 cout << endl << "++ Writing trackLets" << endl;

	trackLets_vector.push_back(m_trackLet);
	m_trackLet.hits.clear();
// cout << " " <<  m_trackLet.trackLetNum << endl;

	    }
 	}
     } 
  }                                

void setTrackLets(trackLet *tr, int hit0_offset, int hit1_offset){
	tr->trackLetNum += 1 ;

	double x0  = hit_Xs[hit0_offset];
        double y0  = hit_Ys[hit0_offset];
        double z0  = hit_Zs[hit0_offset];

	double x1  = hit_Xs[hit1_offset];
        double y1  = hit_Ys[hit1_offset];
        double z1  = hit_Zs[hit1_offset];
	int ze = 800;

	double dx =  (x1 -x0)/(z1 - z0) ;
        double dy =  (y1 -y0)/(z1 - z0) ;

	tr->m_dx =  dx ;
        tr->m_dy =  dy ;
        tr->m_ex = (x1-x0)*(ze-z0)/(z1-z0) + x0;
        tr->m_ey = (y1-y0)*(ze-z0)/(z1-z0) + y0;

	tr->hits.push_back(hit0_offset);
	tr->hits.push_back(hit1_offset);

	//cout << tr->hits[1] << endl;
       }

void trackLetsSelection(vector<trackLet>& vector_trackLets, vector<track>& FinalTracks){
//vector< vector<trackLet> >::iterator iter_ii;
vector<trackLet>::iterator iter_i;
vector<trackLet>::iterator iter_j;
vector<pair <int,int> > compatible1;
vector<pair <int,int> > compatible2;
//clean_vector2(compatible1);
double distdxdy;
double distexey;
for(iter_i=vector_trackLets.begin(); iter_i!=vector_trackLets.end(); iter_i++){

      //  for(iter_j=vector_trackLets.begin(); iter_j!=vector_trackLets.end(); iter_j++){
          for(iter_j=iter_i+1; iter_j!=vector_trackLets.end(); iter_j++){

	//    if(iter_i->trackLetNum != iter_j->trackLetNum){
		distdxdy = pow((iter_i->m_dx - iter_j->m_dx),2)+pow((iter_i->m_dy - iter_j->m_dy),2);
		distexey = pow((iter_i->m_ex - iter_j->m_ex),2)+pow((iter_i->m_ey - iter_j->m_ey),2);
		if(distdxdy <0.00005){
			pair<int,int> comp1(iter_i->trackLetNum,iter_j->trackLetNum); 
			compatible1.push_back(comp1);
//	cout << "TackLet = " << comp1.first << " " << comp1.second << endl;
		}
		if(distexey <50){
                        pair<int,int> comp2(iter_i->trackLetNum,iter_j->trackLetNum);
                        compatible2.push_back(comp2);
        //cout << "TackLet = " << comp2.first << " " << comp2.second << endl;
                }
	   // }
        }
}
//clusters
vector< set<int> > Clusters1;
vector< set<int> > Clusters2;
clean_vector2(Clusters1);
clean_vector2(Clusters2);
set<int> Cand;
set<int> Used;
Used.insert(0);
// findClusterdxdy
for (vector <pair <int,int> >::iterator iter_s = compatible1.begin(); iter_s != compatible1.end(); iter_s++ ){
//	if (value == (*iter_s).first) continue;
	const bool is_in = Used.find((*iter_s).first) != Used.end();
	if (is_in) continue;	
		
       Cand.insert((*iter_s).first);
       Cand.insert((*iter_s).second);
       Used.insert((*iter_s).first);
       Used.insert((*iter_s).second);

     for (vector <pair <int,int> >::iterator iter_v = iter_s +1; iter_v != compatible1.end(); iter_v++){
	if ((*iter_s).first == (*iter_v).first){ 
       Cand.insert((*iter_v).second);
	Used.insert((*iter_v).second);}
//	if((*iter_s).second == (*iter_v).first){
//	Cand.insert((*iter_v).second);
//	Used.insert((*iter_v).second);

//	}
     }
	Clusters1.push_back(Cand);
//for (set<int>::iterator it=Cand.begin(); it!=Cand.end(); ++it){
 //   std::cout << ' ' << *it;}
 // std::cout << '\n';

	Cand.clear();
 }
 Used.clear();
// findClusterexey
        for (vector <pair <int,int> >::iterator iter_s = compatible2.begin(); iter_s != compatible2.end(); iter_s++ ){
        const bool is_in = Used.find((*iter_s).first) != Used.end();
           if (is_in) continue;
       Cand.insert((*iter_s).first);
       Cand.insert((*iter_s).second);
       Used.insert((*iter_s).first);
       Used.insert((*iter_s).second);

     	for (vector <pair <int,int> >::iterator iter_v = iter_s +1; iter_v != compatible2.end(); iter_v++){
        	if ((*iter_s).first == (*iter_v).first){
       		Cand.insert((*iter_v).second);
        	Used.insert((*iter_v).second);}
	}
        Clusters2.push_back(Cand);
	Cand.clear();
 	}
 Used.clear();*/

/*for(vector<set<int> >::iterator iter = Clusters1.begin(); iter != Clusters1.end(); iter++){
	cout << "cluster contains:";
  for (set<int>::iterator it=(*iter).begin(); it!=(*iter).end(); ++it){
    std::cout << ' ' << *it;}
  std::cout << '\n';
}*/

    //do whatever to result
//intersection
/*vector< set<int> > Clusters;
clean_vector2(Clusters);
set<int> v;
//set<int>::iterator it;
	for(vector<set<int> >::iterator iter_c1 = Clusters1.begin(); iter_c1 != Clusters1.end(); iter_c1++){
		for(vector<set<int> >::iterator iter_c2 = Clusters2.begin(); iter_c2 != Clusters2.end(); iter_c2++){
		set_intersection ((*iter_c1).begin(), (*iter_c1).end(), (*iter_c2).begin(), (*iter_c2).end(), inserter(v, v.end()));
                if ( v.size()<= 1 ) continue; 
			Clusters.push_back(v);
			v.clear();
		}
	}
//for(vector<set<int> >::iterator iter = Clusters.begin(); iter != Clusters.end(); iter++){
  //      cout << "cluster contains:" ;
 // for (set<int>::iterator it=(*iter).begin(); it!=(*iter).end(); ++it){
  //  std::cout << ' ' << *it;}
  //std::cout << '\n';
//	}
//convert tracklets to hits
vector<trackLet>::iterator it;
track m_trk;
m_trk.trackHitsNum =0;
//vector< track > FinalTracks;
for(vector<set<int> >::iterator iter = Clusters.begin(); iter != Clusters.end(); iter++){
      for (set<int>::iterator iter_set=(*iter).begin(); iter_set!=(*iter).end(); iter_set++){
	int value = *iter_set;

	for( vector<trackLet>::iterator iter_w=vector_trackLets.begin(); iter_w!=vector_trackLets.end(); iter_w++){
  	if (iter_w->trackLetNum == value){
	
	int hit0 = iter_w->hits[0];
 	int hit1 = iter_w->hits[1];
	setTrack(&m_trk , hit0, hit1);
	//cout<< iter_w->trackLetNum <<endl;
	 }
	}
      }
FinalTracks.push_back(m_trk);
m_trk.trackHitsNum =0;
m_trk.Hits.clear();
  }
//for(vector<track >::iterator iterk = FinalTracks.begin(); iterk != FinalTracks.end(); iterk++){
//	for (set<int>::iterator its=(*iterk).Hits.begin(); its!=(*iterk).Hits.end(); its++){

//cout << " hitSensor " << hit_sensorNums [(*its)];}
//std::cout << '\n';

//}
}//end

void setTrack(track *trk, int hit0, int hit1){
trk->trackHitsNum += 1;
trk->Hits.insert(hit0);
trk->Hits.insert(hit1);
}


template<class T>
void clean_vector(vector<vector<T> >& _tracks_vector){
        for(typename vector<vector<T> >::iterator it = _tracks_vector.begin(); it != _tracks_vector.end(); it++){
                (*it).clear();
        }
}
template<class T>
void clean_vector2(vector<T>& _tracks_vector){
        for(typename vector<T>::iterator it = _tracks_vector.begin(); it != _tracks_vector.end(); it++){
                (*it).clear();
        }
}

template<class T>
void printResultTracks(vector<T> tracks, int event_no, string track_folder_container){
        string track_filename = track_folder_container + "//tracks_" + toString<int>(event_no) + ".txt";
        ofstream track_file(track_filename.c_str());
        track_file << std::setprecision(3);

        int z;

        int t_no = 0;
        for (typename vector<T>::iterator it = tracks.begin(); it != tracks.end(); it++){
 track_file << "track " << t_no++  << std::endl;
               for (set<int>::iterator ith = (*it).Hits.begin(); ith != (*it).Hits.end(); ith++){
                        z = (*ith);
                        track_file << "hit " << hit_ids[(*ith)] << " s " << hit_sensorNums[(*ith)] << " ("
                                           << hit_Xs[(*ith)] << ", " << hit_Ys[(*ith)] << ", " << hit_Zs[(*ith)] << ")" << endl;
                }
                track_file << endl;
        }

        track_file.close();
}*/

                                                        
