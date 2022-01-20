// Implementation of the KMeans Algorithm
// reference: http://mnemstudio.org/clustering-k-means-example-1.htm
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <string>



using namespace std;

class Point
{
private:
	int id_point, id_cluster;
	vector<double> values;
	int total_values;
	string name;

public:
	Point(int id_point, vector<double>& values, string name = "")
	{
		this->id_point = id_point;
		total_values = values.size();

		for(int i = 0; i < total_values; i++)
			this->values.push_back(values[i]);

		this->name = name;
		id_cluster = -1;
	}
	~Point(){}

	int getID()
	{
		return id_point;
	}

	void setCluster(int id_cluster)
	{
		this->id_cluster = id_cluster;
	}

	int getCluster()
	{
		return id_cluster;
	}

	double getValue(int index)
	{
		return values[index];
	}

	int getTotalValues()
	{
		return total_values;
	}

	void addValue(double value)
	{
		values.push_back(value);
	}

	string getName()
	{
		return name;
	}
};

class Cluster
{
private:
	int id_cluster;
	vector<double> central_values;
	vector<Point> points;

public:
	Cluster(int id_cluster, Point point)
	{
		this->id_cluster = id_cluster;

		int total_values = point.getTotalValues();

		for(int i = 0; i < total_values; i++)
			central_values.push_back(point.getValue(i));

		points.push_back(point);
	}
	~Cluster(){
		points.clear();
//		this ->	id_cluster = -1;
//		cout <<"POINTS DELETED  "<<endl;
	}


	void addPoint(Point point)
	{
		points.push_back(point);
	}

	bool removePoint(int id_point)
	{
		int total_points = points.size();

		for(int i = 0; i < total_points; i++)
		{
			if(points[i].getID() == id_point)
			{
				points.erase(points.begin() + i);
				return true;
			}
		}
		return false;
	}

	double getCentralValue(int index)
	{
		return central_values[index];
	}

	void setCentralValue(int index, double value)
	{
		central_values[index] = value;
	}

	Point getPoint(int index)
	{
		return points[index];
	}

	int getTotalPoints()
	{
		return points.size();
	}

	int getID()
	{
		return id_cluster;
	}
};

class KMeans
{
private:
	int K; // number of clusters
	int total_values, total_points, max_iterations;
	vector<Cluster> clusters;

	double getAverageDistanceClusterC(int c)
	{
		int total_points_cluster =  clusters[c].getTotalPoints();
		if(total_points_cluster>0){
			double sum = 0.0;
			double av_dist = 0;
			for(int j = 0; j < total_points_cluster; j++)
			{
				for(int p = 0; p < total_values; p++){
					sum += pow(clusters[c].getCentralValue(p) -
						   clusters[c].getPoint(j).getValue(p), 2.0);
				}
				av_dist += sqrt(sum); // sum SE E' il quadrato della distanza ALTRIMENTI sqrt(sum)!!!
			}
			av_dist = av_dist / total_points_cluster;
			return av_dist;
		}else
		{
			return 0;
		}
	}

	double getAverageDistanceDimPClusterC(int c, int p)
	{
		int total_points_cluster =  clusters[c].getTotalPoints();
		double sum = 0.0;
		double av_dist = 0;
		for(int j = 0; j < total_points_cluster; j++)
		{
			sum += pow(clusters[c].getCentralValue(p) -
			   clusters[c].getPoint(j).getValue(p), 2.0);
			av_dist += sqrt(sum);
		}
		av_dist = av_dist / total_points_cluster;
		return av_dist;
	}





	// return ID of nearest center (uses euclidean distance)
	int getIDNearestCenter(Point point)
	{
		double sum = 0.0, min_dist;
		int id_cluster_center = 0;

		for(int i = 0; i < total_values; i++)
		{
			sum += pow(clusters[0].getCentralValue(i) -
					   point.getValue(i), 2.0);
		}

		min_dist = sqrt(sum);

		for(int i = 1; i < K; i++)
		{
			double dist;
			sum = 0.0;

			for(int j = 0; j < total_values; j++)
			{
				sum += pow(clusters[i].getCentralValue(j) -
						   point.getValue(j), 2.0);
			}

			dist = sqrt(sum);

			if(dist < min_dist)
			{
				min_dist = dist;
				id_cluster_center = i;
			}
		}

		return id_cluster_center;
	}

	vector<double> random_initialization()
	{
		vector<double> clust_in;
		double max_val = 1.0;
		for(int j=0;j<total_values;j++)
		{
			double value = (double)rand() / ((double)RAND_MAX+1.0);
			clust_in.push_back(value);
		}
		return clust_in;
	}






public:
	KMeans(int K, int total_points, int total_values, int max_iterations)
	{
		this->K = K;
		this->total_points = total_points;
		this->total_values = total_values;
		this->max_iterations = max_iterations;
	}
	~KMeans(){
		//cout <<" KMEANS DELETED"<<endl;
		clusters.clear();
	}

	void run(vector<Point> & points, int & count)
	{
		if(count>1){
		clusters.clear();
		}
		for(int i = 0; i < total_points; i++)
		{
			points[i].setCluster(-1);
		}

		if(K > total_points)
			return;

////////////INITIALIZATION WITH RANDOM POINT
		//WRONGGGG
/*		for(int i = 0; i < K; i++) // K = NUMBER OF CLUSTERS
		{
			
			ostringstream convert;
			convert << i;
			string point_name = "Cluster "+convert.str();;
			vector<double> init_vect = random_initialization();
			int int_ind = -1 -1*i;
			Point p(int_ind, init_vect, point_name);
			Cluster cluster(i, p);
			clusters.push_back(cluster);
			clusters[i].removePoint(p.getID());
		}
*/
/*
		for(int i=0;i<K;i++){
			for(int j = 0; j < total_values; j++){
			double value = (double)rand() / ((double)RAND_MAX+1.0);
			clusters[i].setCentralValue(j, value);
			}
		}
*/


////////////////INITIALIZATION CHOOSING CLUSTER POINT FROM EXISTING POINT
		vector<int> prohibited_indexes;
		prohibited_indexes.clear();
		// choose K distinct values for the centers of the clusters
		for(int i = 0; i < K; i++) // K = NUMBER OF CLUSTERS
		{
			while(true)
			{
				int index_point = rand() % total_points; //CHOOSE RANDOM POINT AMONG ALL POINTS
				if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
						index_point) == prohibited_indexes.end())
				{
					prohibited_indexes.push_back(index_point);
					points[index_point].setCluster(i); // SET RANDOM POINT AS INITIAL CENTRAL CLUSTER
					Cluster cluster(i, points[index_point]);
					clusters.push_back(cluster);
					break;
				}
			}
		}

		/* DEBUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUG
		int check2 = 0.0;
		for(int i=0; i<K;i++){
		int total_points_cluster = clusters[i].getTotalPoints();
		cout << "Cluster " << clusters[i].getID();

		check2 += clusters[i].getTotalPoints();
			for(int j = 0; j < total_points_cluster; j++)
			{
				cout << "  Point " << clusters[i].getPoint(j).getID() + 1 << ": ";
				string point_name = clusters[i].getPoint(j).getName();

				if(point_name != "")
					cout << "- " << point_name;

				cout << endl;


			}
		}
		cout <<"CHECK: " << check2 <<endl ;
		*/



////////////////INITIALIZATION CHOOSING CLUSTER POINT FROM EXISTING POINT


		int iter = 1;

		while(true)
		{
			bool done = true;

			// associates each point to the nearest center
			for(int i = 0; i < total_points; i++)
			{	
				//bool aa = true;
				//bool bb = true;
				
				int id_old_cluster = points[i].getCluster();
				int id_nearest_center = getIDNearestCenter(points[i]);
				if(id_old_cluster != id_nearest_center)
				{
					if(id_old_cluster != -1){
					//	cout << "DESTRUCTION OF id " << id_old_cluster << " id point " <<points[i].getID() <<endl;
					//	aa = false;
						clusters[id_old_cluster].removePoint(points[i].getID());
					}
					points[i].setCluster(id_nearest_center);
					clusters[id_nearest_center].addPoint(points[i]);
					//bb = false;
					//cout << "ADDITION OF id " << id_nearest_center << " id point " <<points[i].getID() <<endl;

					done = false;
					//if(aa != bb){cout << "STOP!!!!" << endl;}
				}
			}

//			ofstream output;
//			output.open ("./output/cluster_evolution.csv", ios::out | ios::app); 
	
			// recalculating the center of each cluster
			for(int i = 0; i < K; i++)
			{
				//cout << count << "," << K << "," << i;
//				output << count << "," << K << "," << i;
				for(int j = 0; j < total_values; j++)
				{
					int total_points_cluster = clusters[i].getTotalPoints();
					double sum = 0.0;

					if(total_points_cluster > 0)
					{
						for(int p = 0; p < total_points_cluster; p++)
							sum += clusters[i].getPoint(p).getValue(j);
							clusters[i].setCentralValue(j, sum / total_points_cluster);
				//			cout <<","<< clusters[i].getCentralValue(j);
//							output << ","<<clusters[i].getCentralValue(j);
					}
				}
				//cout <<endl;
//				output<<endl;
			}

//			output.close();	

			if(done == true || iter >= max_iterations)
			{
				/* DEBUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUG
				int check2 = 0.0;
				for(int i=0; i<K;i++){
				check2 += clusters[i].getTotalPoints();
				}
				cout  << "CHECK2: " << check2 <<endl;
				*/
//				cout << "Break in iteration " << iter << "\n\n";
				break;
			}
			iter++;
		}
	}//FINISH FUNCTION RUN





////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////EUCLIDEAN DISTANCE - EVALUATING INDEX/////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////
	double QuantizationError()
	{
		double quanterror = 0.0;
		for(int i = 0; i<K;i++)
		{
			quanterror += getAverageDistanceClusterC(i);
		}
		quanterror = quanterror / K;
		return quanterror;
	}


	double DBIndex_max()
	{
		vector<double> di;
		double sum = 0.0;
		double res = 0.0;
		vector<double> DB;
		for(int i=0;i<K;i++){
			di.push_back(getAverageDistanceClusterC(i));
		}

		for(int i=0;i<K;i++){ //SUM ON ALL CLUSTERS
		for(int j=0;j<K;j++){ //SECOND SUM ON ALL CLUSTERS
		int total_points_clusterA = clusters[i].getTotalPoints();
		int total_points_clusterB = clusters[j].getTotalPoints();

			if(i!=j){  //DA VALUTARE SE INSERIRE && total_points_clusterA > 0 && total_points_clusterB > 0
				sum = 0;
				for(int p = 0; p<total_values;p++){	//SUM ON ALL ATTRIBUTES 
						sum += pow(clusters[i].getCentralValue(p) -
							   clusters[j].getCentralValue(p), 2.0);
				}
				sum = sqrt(sum);
				DB.push_back((di[i]+di[j])/sum);
			}else
			{
				DB.push_back(0);
			}
		}
		res += vect_max(DB);

		}
		return	res / K;
	}



	double SilhouetteIndex2()
	{
		vector<vector <double> > sitot;
		double sisumtot= 0.0;
		for(int k=0;k<K;k++)
		{
		int total_points_clusterA = clusters[k].getTotalPoints();
		vector<double> si;
		double sisum = 0.0;
		if(total_points_clusterA>0)
		{
			for(int a=0;a<total_points_clusterA;a++)
			{
//cout <<"k: "<<k << " a: " << a << endl;
				vector<double> vect_clust;
				for(int kk=0;kk<K;kk++)
				{
					int total_points_clusterB = clusters[kk].getTotalPoints();
					if(k!=kk && total_points_clusterB > 0)
					{
					vect_clust.push_back(ClusterKPointAClusterBCorrelation(a,k,kk));
//cout <<"k: "<<k << " a: " << a <<" kk: "<<kk<< " dzi': " <<ClusterKPointAClusterBCorrelation(a,k,kk)<<endl;
					}
				}
				double first_neigh = vect_min(vect_clust);
				double inn_clust = ClusterKPointACorrelation(a,k);
				double maxd = maxAB(first_neigh,inn_clust);
				si.push_back((first_neigh - inn_clust)/maxd);
				sisum += (first_neigh - inn_clust)/maxd;
//cout <<"d'zi: "<<first_neigh<<" dzi: "<<inn_clust << " sisum: " << sisum << endl;
				vect_clust.clear();
			}
		sisum = sisum / total_points_clusterA;
		sitot.push_back(si);
		sisumtot += sisum;
		}
		si.clear();

		}
		return sisumtot = sisumtot / K;
	}



//RAND INDEX
double Rand_index(vector<Point> & points, vector<Point> & points_real){
	int TPa = 0;
	int TNa = 0;
	double Rand;
	for(int i = 0; i < total_points; i++){
		for(int j=i+1; j < total_points; j++){
			if(points[i].getCluster()==points[j].getCluster() && points_real[i].getCluster()==points_real[j].getCluster())
			{
				TPa++;
			}
			if(points[i].getCluster()!=points[j].getCluster() && points_real[i].getCluster()!=points_real[j].getCluster())
			{
				TNa++;
			}
//			cout <<"i: " << i <<" j: " << j << " " << TPa << " " << TNa <<endl;
		}
	}

	Rand = 2.0*(TPa+TNa) / (total_points*(total_points-1));
//	cout <<"Rand: "<< Rand <<endl;

	return Rand;

}

double Jaccard_index(vector<Point> & points, vector<Point> & points_real){
	int TPa = 0;
	int TNa = 0;
	int FPa = 0;
	int FNa = 0;
	double Jacc = 0.0;
	for(int i = 0; i < total_points; i++){
		for(int j=i+1; j < total_points; j++){
			TPa += TP(i,points,j,points_real);
			TNa += TN(i,points,j,points_real);
			FPa += FP(i,points,j,points_real);
			FNa += FN(i,points,j,points_real);
//		cout <<"i: " << i <<" j: " << j << " " << TPa << " " << TNa << " " <<FPa << " " << FNa <<endl;
		}
	}
	Jacc = (double)TPa/(double)(TPa+TNa+FPa);
//	cout <<"Jacc: "<< Jacc <<endl;
	return Jacc;
}

double Fowlkes_index(vector<Point> & points, vector<Point> & points_real){
	int TPa = 0;
	int FPa = 0;
	int FNa = 0;
	double fowlkes;
	for(int i = 0; i < total_points; i++){
		for(int j=i+1; j < total_points; j++){
			TPa += TP(i,points,j,points_real);
			FPa += FP(i,points,j,points_real);
			FNa += FN(i,points,j,points_real);
//		cout <<"i: " << i <<" j: " << j << " " << TPa  << " " <<FPa << " " << FNa <<endl;
		}
	}
	fowlkes = (double)sqrt((double)TPa / ((double)(TPa+FPa))*(double)TPa/((double)(TPa+FNa))    );
//	cout <<"fow: "<< fowlkes <<endl;
	return fowlkes;	

}

int TP(int i,vector<Point> & points, int j,vector<Point> & points_real){
	if(points[i].getCluster()==points[j].getCluster() && points_real[i].getCluster()==points_real[j].getCluster()){
		return 1;
	}else{
		return 0;
	}
}

int TN(int i,vector<Point> & points,int  j,vector<Point> & points_real){
	if(points[i].getCluster()!=points[j].getCluster() && points_real[i].getCluster()!=points_real[j].getCluster()){
		return 1;
	}else{
		return 0;
	}
}

int FP(int i,vector<Point> & points,int j,vector<Point> & points_real){
	if(points[i].getCluster()==points[j].getCluster() && points_real[i].getCluster()!=points_real[j].getCluster()){
		return 1;
	}else{
		return 0;
	}
}


int FN(int i,vector<Point> & points,int j,vector<Point> & points_real){
	if(points[i].getCluster()!=points[j].getCluster() && points_real[i].getCluster()==points_real[j].getCluster()){
		return 1;
	}else{
		return 0;
	}
}





//FUNCTION TO COMPUTE INTERCORRELATION BETWEEN TWO CLUSTERS
	double ClusterAClusterBCorrelation(int k, int kk)
	{
		int total_points_clusterA = clusters[k].getTotalPoints();
		double sum = 0.0;
		if(total_points_clusterA > 0)
		{
			for(int a=0;a<total_points_clusterA; a++)
			{
				sum += ClusterKPointAClusterBCorrelation(a,k,kk);
			}
		}
		sum = sum / (total_points_clusterA);

		return sum;			
	}

	double ClusterKPointAClusterBCorrelation(int a, int k, int kk)
	{
		int total_points_clusterB = clusters[kk].getTotalPoints();
		double sum = 0.0;
		if(total_points_clusterB > 0)
		{
			for(int b=0;b<total_points_clusterB; b++)
			{
				sum += ClusterKPointAClusterBPointBDistance(a,b,k,kk);
			}
		}
		sum = sum / (total_points_clusterB);
		return sum;				

	}




	double ClusterKPointAClusterBPointBDistance(int a, int b, int k, int kk)
	{
		double sum = 0.0;
		for(int j= 0; j<total_values;j++)
		{
			sum += pow(clusters[k].getPoint(a).getValue(j) -
			clusters[kk].getPoint(b).getValue(j), 2.0);
		}
		return sqrt(sum);
	}



//FUNCTION TO COMPUTE CORRELATION WITHIN A SINGLE CLUSTER
	double ClusterKCorrelation(int k)
	{
		int total_points_cluster = clusters[k].getTotalPoints();
		double sum = 0.0;
		if(total_points_cluster > 0)
		{
			for(int a=0;a<total_points_cluster; a++)
			{
				sum += ClusterKPointACorrelation(a,k);
			}

			if(total_points_cluster == 1)
			{
				return sum / (total_points_cluster);
			}else
			{
				return sum / (total_points_cluster-1);
			}
		}else
		{
			return 0;
		}
	}

	double ClusterKPointACorrelation(int a, int k)
	{
		double sum = 0.0;
		int total_points_cluster = clusters[k].getTotalPoints();
		for(int b=0; b < total_points_cluster;b++)
		{
				sum += ClusterKPointAPointBDistance(a,b,k);
		}
		if(total_points_cluster == 1)
		{
			return sum / (total_points_cluster);
		}else{
			return sum / (total_points_cluster-1);
		} 
	}


	double ClusterKPointAPointBDistance(int a, int b, int k)
	{
		double sum = 0.0;
		for(int j= 0; j<total_values;j++)
		{
			sum += pow(clusters[k].getPoint(a).getValue(j) -
			clusters[k].getPoint(b).getValue(j), 2.0);
		}
		return sqrt(sum);
	}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////PRINT FUNCTION////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////


	//PRINT EVALUATING INDICES
	void print_indices(vector<Point> & points, int & count)
	{
		ofstream output;
		output.open ("./output/indices.csv", ios::out | ios::app); 
		output << count << " " << K << " " 
		<< QuantizationError() 
		<<" " << DBIndex_max()
		<<" "<<SilhouetteIndex2()<< endl;
		output.close();
	}



	void print(vector<Point> & points, int & count){
		// shows elements of clusters on terminal
		/*
		for(int i = 0; i < K; i++)
		{
			int total_points_cluster =  clusters[i].getTotalPoints();
					cout << "Cluster " << clusters[i].getID() + 1 << endl;
			for(int j = 0; j < total_points_cluster; j++)
			{
				cout << "Point " << clusters[i].getPoint(j).getID() + 1 << ": ";
				for(int p = 0; p < total_values; p++)
					cout << clusters[i].getPoint(j).getValue(p) << " ";

				string point_name = clusters[i].getPoint(j).getName();

				if(point_name != "")
					cout << "- " << point_name;

				cout << endl;
			}

			cout << "Cluster values: ";

			for(int j = 0; j < total_values; j++)
				cout << clusters[i].getCentralValue(j) << " ";

			cout << "\n\n";
		}
		*/
		//DEBUG
		int check_number = 0;
		for(int i=0;i<K;i++){
			check_number += clusters[i].getTotalPoints();
		}
		if(check_number!= total_points){
			cout <<endl<<endl<<"WRONG NUMBER OF POINTS "
			<< "right number: " << total_points
			<< "wrong number: " << check_number
			<< "K: "<< K << "count " << count<<endl<<endl;
		}

		ofstream output;
		output.open ("./output/cluster.csv", ios::out | ios::app);
		output << ">>>>>>>>>>>>>>>>>>    "
		<<"Simulazione: " << count << " K: "<< K
		<< "    <<<<<<<<<<<<<<<<<" <<endl << "CHECK: " << check_number
		<< endl << endl;

		for(int i = 0; i < K; i++)
		{
			int total_points_cluster =  clusters[i].getTotalPoints();
			output << "Cluster " << clusters[i].getID() + 1 << endl;
			for(int j = 0; j < total_points_cluster; j++)
			{
				output << "Point " << clusters[i].getPoint(j).getID() + 1 << ": ";
				for(int p = 0; p < total_values; p++)
					output << clusters[i].getPoint(j).getValue(p) << " ";

				string point_name = clusters[i].getPoint(j).getName();

				if(point_name != "")
					output << "- " << point_name;

				output << endl;
			}

			output << "Cluster values: ";

			for(int j = 0; j < total_values; j++)
				output << clusters[i].getCentralValue(j) << " ";

			output << "\n\n";
		}
		output.close();
	}//FINISH PRINT


	void print_terminal(vector<Point> & points){
		for(int i = 0; i < K; i++)
		{
			cout << "Cluster "<<i<<" values: ";
			for(int j=0; j<total_values;j++){
				cout << clusters[i].getCentralValue(j)<<" ";
				
			}
			cout << endl;
		}
	}//FINISH FUNCTION PRINT2


	void print_output(vector<Point> & points, int & count){
		// shows elements of clusters on terminal
		ofstream output;
		output.open ("./output/cluster_tot.csv", ios::out | ios::app);
		for(int i = 0; i < K; i++)
		{
			int total_points_cluster =  clusters[i].getTotalPoints();

			for(int j = 0; j < total_points_cluster; j++)
			{
				output << count << "," << K
				<<","<< clusters[i].getID() + 1
				<<","<< clusters[i].getPoint(j).getID() + 1
				<<","<< clusters[i].getPoint(j).getName();

				for(int p = 0; p < total_values; p++){
					output << ","<<clusters[i].getPoint(j).getValue(p) ;
				}
				output << endl;
			}
		}
		output.close();


		ofstream output_cumulative;
		output_cumulative.open("./output/cluster_cumulative.csv", ios::out | ios::app);
		for(int i = 0; i < K; i++){
			output_cumulative << count << "," << K
			<<","<< clusters[i].getID() + 1;

			for(int j = 0; j < total_values; j++)
				output_cumulative << ","<<clusters[i].getCentralValue(j);
				output_cumulative << endl;
		}
		output_cumulative.close();
	}//FINISH PRINT 3


	void print_finalcluster(vector<Point> & points){
		// shows elements of clusters on terminal
		ofstream output;
		output.open ("./output/final_cluster.csv");
		output << "Id Cluster,Id Point, Point Name";
		for(int i=0;i<total_values;i++){
			output <<",attribute"<<i;
		}
		output<<endl;
		output.close();

		output.open ("./output/final_cluster.csv", ios::out | ios::app);
		for(int i = 0; i < K; i++)
		{
			int total_points_cluster =  clusters[i].getTotalPoints();

			for(int j = 0; j < total_points_cluster; j++)
			{
				output << clusters[i].getID() + 1
				<<","<< clusters[i].getPoint(j).getID() + 1
				<<","<< clusters[i].getPoint(j).getName();

				for(int p = 0; p < total_values; p++){
					output << ","<<clusters[i].getPoint(j).getValue(p) ;
				}
				output << endl;
			}
		}
	}//FINISH PRINT 4


	void print_finalcentroid(vector<Point> & points){
		ofstream output;

		output.open("./output/final_centroid.csv");
		output << "ID_Cluster";
		for(int i=0;i<total_values;i++){
			output <<",attribute"<< i << ",attr min" << i << ",attr max" << i;
		}
		output<<endl;
		output.close();

		output.open("./output/final_centroid.csv", ios::out | ios::app);
		for(int i = 0; i < K; i++){
			output << clusters[i].getID() + 1;
			for(int j = 0; j < total_values; j++){
				double centr_value = clusters[i].getCentralValue(j);
				double variance = getAverageDistanceDimPClusterC(i,j);
				output << ","<< centr_value << "," 
				<< centr_value - variance << ","
				<< centr_value + variance;
			}
				output << endl;
		}
		output.close();
	}


	void print_finalcentroid2(vector<Point> & points, vector<Point> & points_norm){

		for(int i = 0; i < total_points; i++)
		{	
			points[i].setCluster(points_norm[i].getCluster());
		}



		ofstream output;

		output.open("./output/final_centroid2.csv");
		output << "ID_Cluster";
		for(int i=0;i<total_values;i++){
			output <<",attribute"<< i << ",attr min" << i << ",attr max" << i;
		}
		output<<endl;
		output.close();

		output.open("./output/final_centroid2.csv", ios::out | ios::app);
		for(int i = 0; i < K; i++){
			output << clusters[i].getID() + 1;
			for(int j = 0; j < total_values; j++){
				double centr_value = clusters[i].getCentralValue(j);
				double variance = getAverageDistanceDimPClusterC(i,j);
				output << ","<< centr_value << "," 
				<< centr_value - variance << ","
				<< centr_value + variance;
			}
				output << endl;
		}
		output.close();
	}








	void print_finalindices(vector<Point> & points)
	{
		ofstream output;
		output.open ("./output/final_indices.csv"); 
		output << "QuantError"<< " "<<"DBIndex" <<" "<<"Silhouette"<< endl;
		output.close();
		output.open ("./output/final_indices.csv", ios::out | ios::app); 
		output << QuantizationError() 
		<<" " << DBIndex_max()
		<<" "<<SilhouetteIndex2()<< endl;
		output.close();
	}

	void print_final_externalindices(vector<Point> & points,vector <Point> & points_real)
	{
		ofstream output;
		output.open ("./output/final_external_indices.csv"); 
		output << "Rand"<< " "<<"Jacc" <<" "<<"Fowlkes"<< endl;
		output.close();
		output.open ("./output/final_external_indices.csv", ios::out | ios::app); 
		output 
		      << Rand_index(points,points_real) 
		<<" " << Jaccard_index(points,points_real)
		<<" " << Fowlkes_index(points,points_real)
		<< endl;
		output.close();
	}




	double vect_min(vector <double> v){
	    double mic = v[0];
	    for(int i = 0; i != v.size(); ++i)
	    {
	        if(v[i] < mic)
	        mic = v[i];
	    }
	    return mic;
	}

	double vect_max(vector <double> v){
	    double mic = v[0];
	    for(int i = 0; i != v.size(); ++i)
	    {
	        if(v[i] > mic)
	        mic = v[i];
	    }
	    return mic;
	}

	double maxAB(double a, double b){
		if(a>b){
			return a;	
		}else{
			return b;
		}
	}

	double minAB(double a, double b){
		if(a<b){
			return a;	
		}else{
			return b;
		}
	}
	//EUCLIDEAN DISTANCE between two point 1Dim
	double ABdist(double a, double b)
	{
		if(a>=b){
			return a-b;
		}
		else{
			return b-a;
		}
	}

	double ABvectDist(vector<double> a,vector<double> b){
			double sum = 0.0;
			double av_dist = 0;
			for(int i = 0; i < a.size(); i++)
			{
				sum += pow(a[i] - b[i], 2.0);
				av_dist += sqrt(sum);
			}
			av_dist = sqrt(sum);
			return av_dist;
	}

}; //END OF CLASS KMEANS

class Data_temp {
   public:
      vector <vector <string> > temp; 
};


double vect_min(vector <double> v){
    double mic = v[0];
    for(int i = 0; i != v.size(); ++i)
    {
        if(v[i] < mic)
        mic = v[i];
    }
    return mic;
}

double vect_max(vector <double> v){
    double mic = v[0];
    for(int i = 0; i != v.size(); ++i)
    {
        if(v[i] > mic)
        mic = v[i];
    }
    return mic;
}

double maxAB(double a, double b){
	if(a>b){
		return a;	
	}else{
		return b;
	}
}

double minAB(double a, double b){
	if(a<b){
		return a;	
	}else{
		return b;
	}
}
//EUCLIDEAN DISTANCE between two point 1Dim
double ABdist(double a, double b)
{
	if(a>=b){
		return a-b;
	}
	else{
		return b-a;
	}
}

double ABvectDist(vector<double> a,vector<double> b){
		double sum = 0.0;
		double av_dist = 0;
		for(int i = 0; i < a.size(); i++)
		{
			sum += pow(a[i] - b[i], 2.0);
			av_dist += sqrt(sum);
		}
		av_dist = sqrt(sum);
		return av_dist;
	}





int main(int argc, char *argv[])
{

//READ THE FILE WITH THE DATA
	Data_temp results1;

	ifstream infile2( "./datasets/results2.csv" ); //2 DATA
	while (infile2){
	    string s2;
	    if (!getline( infile2, s2 )) break;

	    istringstream ss2( s2 );
	    vector <string> record_testo;

	    while (ss2)
	    {
	      string s2;
	      if (!getline( ss2, s2, ',' )) break;
	      record_testo.push_back( s2 );
	    }
   	    results1.temp.push_back( record_testo );

	  }
	  if (!infile2.eof())
	  {
	    cerr << "Fooey!\n";
	  }
//DEBUG AND CHECK
/*
	for (int i = 0; i < results1.temp.size(); i++) {
		for (int j = 0 ; j<  results1.temp[0].size(); j++){
		    cout <<"i: "<< i << " j: " << j <<" " << results1.temp[i][j] << "\n";
		}
	}
*/




//READ DATA FOR CLUSTER
	Data_temp id_temp;

	ifstream infile3( "./datasets/id_cluster.csv" );
	while (infile3){
	    string s2;
	    if (!getline( infile3, s2 )) break;

	    istringstream ss2( s2 );
	    vector <string> record_testo;

	    while (ss2)
	    {
	      string s2;
	      if (!getline( ss2, s2, ',' )) break;
	      record_testo.push_back( s2 );
	    }
   	    id_temp.temp.push_back( record_testo );

	  }
	  if (!infile3.eof())
	  {
	    cerr << "Fooey!\n";
	  }





//SET DATA FOR K-MEANS
	srand (time(NULL));
	int total_points, total_values, K, max_iterations, has_name;
	// total_points -> TOT OF DATA
	// total_values -> TOT OF ATTRIBUTES
	// K -> N° of CLUSTERS
	// max_iterations -> stops by max numb of iterations or if no data point exchange cluster.
	// has_name -> label of point

	total_points = results1.temp.size();
	total_values = results1.temp[0].size()-1;
	K = 4;
	max_iterations = 1000;
	vector<Point> points;
	string point_name;


	cout << "n° of data: " << total_points<<endl;
	cout << "n° of attributes: "<< total_values<<endl<<endl;


//FIRST RUN
	for(int i = 0; i < total_points; i++)
	{
		vector<double> values;
		point_name = results1.temp[i][0];
		cout << point_name<<" ";
		for(int j = 1; j <= total_values; j++)
		{
				double value;
				value = atof(results1.temp[i][j].c_str());
				cout << value<<" ";
				values.push_back(value);
		}
		cout <<endl;
		Point p(i, values, point_name);
		points.push_back(p);
	}
//FIRST RUN TO SET ID CLUSTER
	vector<Point> points_id;

	for(int i = 0; i < total_points; i++)
	{
		vector<double> values;
		point_name = id_temp.temp[i][0];
		cout << point_name<<" ";
		for(int j = 1; j <= total_values; j++)
		{
				double value = 0;
				cout << value<<" ";
				values.push_back(value);
		}
		cout <<endl;
		Point p(i, values, point_name);
		points_id.push_back(p);
		points_id[i].setCluster(atof(results1.temp[i][2].c_str()));
	}



//SET NORMALIZE POINTS
	vector<Point> norm_points;
	vector<double> val_min,val_max;
 	for(int i = 1; i <= total_values; i++)
	{
		vector<double> values;
		for(int j = 0; j < total_points; j++)
		{
				double value;
				value = atof(results1.temp[j][i].c_str());
				values.push_back(value);
		}
		val_min.push_back(vect_min(values));
		val_max.push_back(vect_max(values));
	}

	ofstream output_norm;
	output_norm.open ("./output/data_norm.csv");
	output_norm << "Name,kWh/m2,day/night,kWh" << endl;
	for(int i = 0; i < total_points; i++)
	{
		vector<double> values;
		point_name = results1.temp[i][0];
		output_norm << point_name;
		for(int j = 1; j <= total_values; j++)
		{
				double value;
				value = (atof(results1.temp[i][j].c_str())-val_min[j-1])/(val_max[j-1]-val_min[j-1]);
				output_norm <<","<< value;
				values.push_back(value);
		}
		output_norm <<endl;
		Point p(i, values, point_name);
		norm_points.push_back(p);
	}
   output_norm.close();





//MONTECARLO 1000 simulations! 
/*
	KMeans kmeans(K, total_points, total_values, max_iterations);
	int a = 0;
	kmeans.run(points,a);
	kmeans.print(points,a);
//PRINT FIRST LINE IN OUTPUT FILE AND ERASE PREVIOUS CONTENT
	ofstream output;
	output.open ("./output/cluster_tot.csv");
	output << "Iteration,Id Cluster,Id Point, Point Name";
	for(int i=0;i<total_values;i++){
		output <<",attribute"<<i;
	}
	output<<endl;
	output.close();

	ofstream output_cumulative;
	output_cumulative.open("./output/cluster_cumulative.csv");
	output_cumulative << "Iteration,Id Cluster";
	for(int i=0;i<total_values;i++){
		output_cumulative <<",attribute"<<i;
	}
	output_cumulative<<endl;
	output_cumulative.close();
	for(int i = 1; i < 1000; i++){
		cout << i << endl;
		kmeans.run(points);
		kmeans.print_terminal(points);
		kmeans.print_output(points, i);
	}
*/




//PRINT FIRST LINE FOR EACH FILE OF CLUSTER
	ofstream output;
	output.open ("./output/indices.csv"); 
	output <<"Iteration"<<" "<< "K" << " " << "QuantError"<< " "<<"DBIndex" <<" "<<"Silhouette"<< endl;
	output.close();


	output.open ("./output/cluster_tot.csv");
	output << "Iteration,K,Id Cluster,Id Point, Point Name";
	for(int i=0;i<total_values;i++){
		output <<",attribute"<<i;
	}
	output<<endl;
	output.close();


	output.open("./output/cluster_cumulative.csv");
	output << "Iteration,K,Id Cluster";
	for(int i=0;i<total_values;i++){
		output <<",attribute"<<i;
	}
	output<<endl;
	output.close();


	output.open ("./output/cluster.csv"); 
	output << "ALL DATA" << endl <<endl;
	output.close();


	output.open("./output/cluster_evolution.csv");
	output << "Iteration,K,ID Cluster,";
	for(int i=0;i<total_values;i++){
		output <<",attribute"<<i;
	}
	output<<endl;
	output.close();
	
	output.open ("./output/knee.csv");
	double Qval = 10;

//	vector <vector <vector <double> > >final_indices3;
	for(int n = 9; n < 10; n++){
//		vector<vector <double> > final_indices2;
		vector <double> Qvect,DBvect,Silvect,Rand_vect,Jacc_vect,fow_vect;
		double DBval,Silval,Randval,Jaccval,fowval;
		Jaccval = 0;
		Randval = 0;
		fowval = 0;
		DBval = 10;
		Silval =0;
		Qval = 10;
		for(int i = 1; i < 10; i++){

			K = n;
//			cout << "SIMULAZIONE " << i << " K: " << K << endl;
			KMeans kmeans(K, total_points, total_values, max_iterations);

			kmeans.run(norm_points,i);
//			kmeans.print_terminal(norm_points);
//			kmeans.print_indices(norm_points,i);
//			kmeans.print_output(norm_points, i);
//			kmeans.print(norm_points,i);

//INTERNAL INDICES
			Qvect.push_back(kmeans.QuantizationError());
			DBvect.push_back(kmeans.DBIndex_max());
			Silvect.push_back(kmeans.SilhouetteIndex2());
//EXTERNAL INDICES
			Rand_vect.push_back(kmeans.Rand_index(norm_points,points_id));
			Jacc_vect.push_back(kmeans.Jaccard_index(norm_points,points_id));
			fow_vect.push_back(kmeans.Fowlkes_index(norm_points,points_id));

			if(Qval > kmeans.QuantizationError()){
//			if(Jaccval < kmeans.Jaccard_index(norm_points,points_id)){
				//INTERNAL INDEX
				Qval = kmeans.QuantizationError();
				DBval = kmeans.DBIndex_max();
				Silval = kmeans.SilhouetteIndex2();
				//EXTERNAL INDEX
				Randval = kmeans.Rand_index(norm_points,points_id);
				Jaccval = kmeans.Jaccard_index(norm_points,points_id);
				fowval = kmeans.Fowlkes_index(norm_points,points_id);
			}

//			vector <double> final_indices;
//			final_indices.push_back(kmeans.SilhouetteIndex2());
//			final_indices.push_back(kmeans.DBIndex_max());
//			final_indices.push_back(kmeans.QuantizationError());
//			final_indices2.push_back(final_indices);
//			final_indices.clear();
		}
///////////////////////////////////////////////////////////////////////////////////////////////////////
		//RIFAAREEEEEEEEEEEEEE QUESTO PEZZOOOOOOOOOOOOOOOOOOOOOOOOO///////////////////////
		//SONO SCORRELATIIIIIIIIIII COSI!!!!!!!!!!!!!!!!!!!///////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

		cout 
		<< "FIRST METHOD "<< endl
		<< " K: " << K
		<< " QE: " << vect_min(Qvect) 
		<< " DB: " << vect_min(DBvect) 
		<< " Sil: "<< vect_max(Silvect)
		<< " Rand: " << vect_max(Rand_vect)
		<< " Jacc: " << vect_max(Jacc_vect)
		<< " Fow: " << vect_max(fow_vect)
		<< endl;

		Qvect.clear();
		DBvect.clear();
		Silvect.clear();
		Rand_vect.clear();
		Jacc_vect.clear();

		output.open ("./output/knee.csv", ios::out | ios::app);
		cout 
		<< "SECOND METHOD " <<endl
		<< " K: " << K
		<< " QE: " << Qval
		<< " DB: " << DBval
		<< " Sil: "<< Silval
		<< " Rand: " << Randval
		<< " Jacc: " << Jaccval
		<< " Fowl: " << fowval
		<< endl;
		output
		<< K << " " << Qval << " " << DBval << " " << Silval << " " << Randval << " " << Jaccval << endl;
		output.close();

//		final_indices3.push_back(final_indices2);
//		final_indices2.clear();
	}



K=9;
	double Qval_temp = 10;
	double Jacc_temp = 0;
	for(int i = 1; i < 1000; i++){
		KMeans kmeans(K, total_points, total_values, max_iterations);
		kmeans.run(norm_points,i);
	
//	if(Qval_temp > kmeans.QuantizationError()){
	if(Jacc_temp < kmeans.Jaccard_index(norm_points,points_id)){
//		Qval_temp = kmeans.QuantizationError();
		Jacc_temp = kmeans.Jaccard_index(norm_points,points_id);
		kmeans.print_finalcluster(norm_points);
		kmeans.print_finalcentroid(norm_points);	
		kmeans.print_finalindices(norm_points);
		kmeans.print_final_externalindices(norm_points,points_id);
//		kmeans.print_finalcentroid2(points,norm_points);
	}
	}

//RAND FUNCTION 









	return 0;
}
