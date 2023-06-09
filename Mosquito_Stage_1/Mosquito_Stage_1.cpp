/*############################################################
notes and important points






*/





#include<iostream>
#include<vector>
#include<string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>
#define SEED_NONE 0 
#include <chrono> // Precise time reporting
#include<vector>
#include <algorithm>
#include <sstream>
#include <fstream>
//#include "func.cpp"//



using namespace std;
using std::cout;

//hello

int mosq_count = 0;
int pause_for_interupt;
int file_counter = 0;
int global_zero_mosq=0;

struct mosq {
	int mosq_id;
	short int status; //stage of mosquito like larve, eggs, male... egg=0, larve=1, pupa=2, adult=3
	bool is_female; //is mosquito male or female
	bool is_treated; //if the mosquito is treated or not(wild or treated)
	bool is_mated;	//if the mosquito is ever mated
	short int gene_stat; //checks status of genes of mosquito like DD(Treated)=2, DR(gen_1)=1, RR(Wild)=0
	struct mosq* next;	//pointer for next location
	//struct mosq* prev;  //for double linked list
};

vector <mosq*> male_vector;
vector <mosq*> egg_vector;
vector <mosq*> larve_vector;
vector <mosq*> pupe_vector;
vector <mosq*> female_vector;

vector <int> male_count;
vector <int> female_count;
vector <int> egg_count;
vector <int> larve_count;
vector <int> pupe_count;
vector <int> time_count;

vector <int> DD_mate_count_vect;
vector <int> DR_mate_count_vect;
vector <int> RR_mate_count_vect;


vector <int> larve_ids;
vector <int> egg_ids;

vector <int> Male_DD;
vector <int> Male_count_DD;
vector <int> Male_DR;
vector <int> Male_count_DR;
vector <int> Male_RR;
vector <int> Male_count_RR;

vector <int> eggs_female_vector;
vector <int> eggs_male_RR_vector;
vector <int> eggs_male_DR_vector;

vector <double> female_ratio;
vector <double> male_ratio;
vector <double> larve_ratio;
vector <double> pupe_ratio;



//mosq* curr_mosq;
struct mosq* prev_mosq;
//curr_mosq = head;


struct mosq* curr_mosq;
struct mosq* last_mosq;

struct mosq* head = NULL; //the head of the mosquito list
struct mosq* temp = NULL;	//the tail of mosquito list
struct mosq* next_mosq = NULL;

struct mosq_vars {
	int initial_mosqs = 40;
	int total_time = 2; //initial value to loop over time
	unsigned int k = 1;	// randomness of initial mosquito death
	unsigned int mate_oc, next_stg_oc;
	unsigned int treated_oc, gender_oc, death_oc, find_male_oc, linear_egg_oc;
	unsigned long int male_rndm_mosq, no_of_male, max_larve_cap = 15; //caring capacity larve cap 10000
	float mate_prob = 0.9;
	float treated_prob = 0.8, female_death_prob = 0.2, male_death_prob = 0.2, egg_death_prob = 0.2, larve_death_prob = 0.2, pupe_death_prob = 0.2, initial_mate_prob = 0.001;
	float gender_prob = 0.5, death_prob = 0.2, find_male_prob = 0.8, next_stg_prob = 0.9, egg_next_stg_prob = 0.9, larve_next_stage_prob = 0.9, pupe_next_stage_prob = 0.9;
	int death_count = 0;
	int intrvl_of_days_to_add_treated_mosq,num_of_days_to_expmnt, num_of_times_in_each_day;
	int num_of_mosq_to_add;
	double step_size;
	int equlibrium_days, max_egg_cap, percent_males_DD_add, adult_males_on_equibrilium = 0, DR_mate_cntr = 0, DD_mate_cntr = 0, RR_mate_cntr = 0, linear_egg_cap=0,global_egg_count=0;
	//input the data from file
	//death prob for each status
	//nxt stg fore each status
};

mosq_vars env_1;


void read_variables() //function for reading varible values from CSV files
{
	string fname = "variables.csv";
	vector<string> vals;
	vector<string> vals1;
	string line, word;

	fstream file(fname, ios::in);
	if (file.is_open())
	{
		while (getline(file, line))
		{
			stringstream str(line);
			while (getline(str, word, ','))
				vals.push_back(word);
		}
	}
	else
		cout << "Could not open the file\n";

	string k1 = vals.at(1);
	int k2 = stoi(k1);
	cout << "\n k1 is : " << k2;
	cout << k2 + 3;



	//assigning values to env variables

	env_1.initial_mosqs = stoi(vals.at(1));
	env_1.max_larve_cap = stoi(vals.at(3));
	env_1.mate_prob = stof(vals.at(5));
	env_1.treated_prob = stof(vals.at(7));
	env_1.female_death_prob = stof(vals.at(9));
	env_1.male_death_prob = stof(vals.at(11));
	env_1.egg_death_prob = stof(vals.at(13));
	env_1.larve_death_prob = stof(vals.at(15));
	env_1.pupe_death_prob = stof(vals.at(17));
	env_1.gender_prob = stof(vals.at(19));
	env_1.find_male_prob = stof(vals.at(21));
	env_1.egg_next_stg_prob = stof(vals.at(23));
	env_1.larve_next_stage_prob = stof(vals.at(25));
	env_1.pupe_next_stage_prob = stof(vals.at(27));
	env_1.intrvl_of_days_to_add_treated_mosq = stoi(vals.at(29));
	env_1.percent_males_DD_add = stoi(vals.at(31));
	env_1.num_of_days_to_expmnt = stoi(vals.at(33));
	env_1.num_of_times_in_each_day = stoi(vals.at(35));
	env_1.equlibrium_days = stoi(vals.at(37));
	env_1.max_egg_cap = stoi(vals.at(39));
	env_1.linear_egg_cap = stoi(vals.at(41));




	env_1.step_size = 0.1;// to calculate step size


}

void add_mosq(int new_status, bool is_female_new, bool is_treated_new, bool is_mated_new, int gene_stat_new) { //adds mosquioes or eggs to the environment
	mosq* new_mosq = new mosq;
	mosq_count += 1;
	new_mosq->mosq_id = mosq_count;
	new_mosq->status = new_status;
	new_mosq->is_female = is_female_new;
	new_mosq->is_treated = is_treated_new;
	new_mosq->is_mated = is_mated_new;
	new_mosq->gene_stat = gene_stat_new;

	new_mosq->next = NULL;

	//updating each vector of their own status
	switch (new_status) {
	case 0:
		egg_vector.push_back(new_mosq);
		break;
	case 1:
		larve_vector.push_back(new_mosq);
		break;
	case 2:
		pupe_vector.push_back(new_mosq);
		break;
	case 3:
		if (is_female_new == true) {
			female_vector.push_back(new_mosq);
		}
		else {
			male_vector.push_back(new_mosq);
		}
		break;
	}

	if (head == NULL) { //adding mosquito to the head if it is first mosq
		head = new_mosq;
		last_mosq = new_mosq;
		//last_mosq = head;
	}
	else { //appending mosquitoes
		last_mosq->next = new_mosq;
		last_mosq = new_mosq;

	}

	//new_mosq->prev = head;  // for double linked list
	//head = new_mosq;
}



void prnt_vectr(vector <mosq*> inp_vect, string gender) { //printing vectors fo each type of mosquitoes
	/*cout << "\n" << gender<<" vector size = " << inp_vect.size()<<"\n";
	for (int i = 0; i < inp_vect.size(); i++) {
		cout << inp_vect.at(i)->mosq_id << " ";
	}*/

}

void prnt_all_vect() { // printing all types at once in vector
	/*prnt_vectr(male_vector, "Male");
	prnt_vectr(female_vector, "Female");
	prnt_vectr(egg_vector, "Egg");
	prnt_vectr(larve_vector, "Larve");
	prnt_vectr(pupe_vector, "Pupe");*/
}


void show_mosqs() { //printing mosquito linked list
	//struct mosq* ptr;
	//cout << "\n ********** entered mosquito showwww ***********";
	//ptr = head;
	//int total_pop = egg_vector.size() + larve_vector.size() + pupe_vector.size() + male_vector.size() + female_vector.size(); //finding the population of all mosquitoes
	//cout << "\n the total population of all mosquitoes is: " << total_pop;
	//cout << "\n the number of male mosquitos are: " << male_vector.size();
	//cout << "\n the number of female mosquitos are: " << female_vector.size();
	//cout << "\n the number of egg mosquitos are: " << egg_vector.size();
	//cout << "\n the number of lavre mosquitos are: " << larve_vector.size();
	//cout << "\n the number of pupe mosquitos are: " << pupe_vector.size();
	//if (head == NULL) {
	//	cout << "\n linked list is empty";
	//	return;
	//}

	//cout << "\n the mosquitos added so --------------------------- far are \n";
	//cout << "\n Mosquito ID: ";
	//cout << " Mosquito stage: ";
	//cout << " is_female ";
	//cout << " treated status: ";
	//cout << " mated status: ";
	//while (ptr != NULL) {
	//	cout << "\n " << ptr->mosq_id;
	//	cout << "\t\t " << ptr->status;
	//	cout << "\t\t " << ptr->is_female;
	//	cout << "\t\t " << ptr->is_treated;
	//	cout << "\t\t " << ptr->is_mated;
	//	ptr = ptr->next;
	//}

	//cout << "\n " << head->mosq_id;
	//cout << "\n " << last_mosq->mosq_id;
	//cout << "\n ********** mosquito show completed ***********";
}

std::string doubleToString(double value) {
	std::ostringstream ss;
	ss << value;
	return ss.str();
}

void write_double_CSV(const std::vector<double>& data, const std::string& filename) {
	std::ofstream file(filename);

	// Write the header
	file << "Values" << std::endl;

	// Write the data
	for (const auto& value : data) {
		file << doubleToString(value) << std::endl;
	}

	file.close();
}


void write_csv(std::string filename, std::vector<std::pair<std::string, std::vector<int>>> dataset) {
	// Make a CSV file with one or more columns of integer values
	// Each column of data is represented by the pair <column name, column data>
	//   as std::pair<std::string, std::vector<int>>
	// The dataset is represented as a vector of these columns
	// Note that all columns should be the same size

	// Create an output filestream object
	std::ofstream myFile(filename);

	// Send column names to the stream
	for (int j = 0; j < dataset.size(); ++j)
	{
		myFile << dataset.at(j).first;
		if (j != dataset.size() - 1) myFile << ","; // No comma at end of line
	}
	myFile << "\n";

	// Send data to the stream
	for (int i = 0; i < dataset.at(0).second.size(); ++i)
	{
		for (int j = 0; j < dataset.size(); ++j)
		{
			myFile << dataset.at(j).second.at(i);
			if (j != dataset.size() - 1) myFile << ","; // No comma at end of line
		}
		myFile << "\n";
	}

	// Close the file
	myFile.close();
}

void write_d_csv(std::string filename, std::vector<std::pair<std::string, std::vector<double>>> dataset) {
	// Make a CSV file with one or more columns of integer values
	// Each column of data is represented by the pair <column name, column data>
	//   as std::pair<std::string, std::vector<int>>
	// The dataset is represented as a vector of these columns
	// Note that all columns should be the same size

	// Create an output filestream object
	std::ofstream myFile(filename);

	// Send column names to the stream
	for (int j = 0; j < dataset.size(); ++j)
	{
		myFile << dataset.at(j).first;
		if (j != dataset.size() - 1) myFile << ","; // No comma at end of line
	}
	myFile << "\n";

	// Send data to the stream
	for (int i = 0; i < dataset.at(0).second.size(); ++i)
	{
		for (int j = 0; j < dataset.size(); ++j)
		{
			myFile << dataset.at(j).second.at(i);
			if (j != dataset.size() - 1) myFile << ","; // No comma at end of line
		}
		myFile << "\n";
	}

	// Close the file
	myFile.close();
}




static gsl_rng* _RNG_P = NULL;


//To allocate and initialize:
int alloc_and_setup_gsl_random_generator() {
	gsl_rng_env_setup(); // Gathers environment variables.
	_RNG_P = gsl_rng_alloc(gsl_rng_default); // Allocates the rnd generator.
	if (!_RNG_P) { return 1; }
	return 0;
}

using int_gsl_seed = unsigned long int;


uint64_t get_time_microsseconds()
{
	uint64_t us = std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::high_resolution_clock::now().time_since_epoch())
		.count();
	return us;

}

int_gsl_seed seed_gsl_random_generator(int_gsl_seed seed) {

	if (seed == SEED_NONE) {
		seed = (int_gsl_seed)get_time_microsseconds();
		// seed = time(0); // In seconds; not enough time resolution.
	}
	gsl_rng_set(_RNG_P, seed);
	return seed;

}


void will_mosq_die() {
	float status_death_prob;
	int m_die_vect_cntr = -1;
	int f_die_vect_cntr = -1;
	int e_die_vect_cntr = -1;
	int l_die_vect_cntr = -1;
	int p_die_vect_cntr = -1;
	int death_cntr = 0;
	int total_pop;

		while (temp != NULL) {
			total_pop = egg_vector.size() + larve_vector.size() + pupe_vector.size() + male_vector.size() + female_vector.size();
			if (total_pop < 2) {
				global_zero_mosq = 1;
				cout << "\n\n\n\n zero mosqutoes alive in 1st try: \n\n\n\n";
				return;
			}

			switch (curr_mosq->status) {
			case 0:
				e_die_vect_cntr++;
				status_death_prob = env_1.egg_death_prob * env_1.step_size;
				break;
			case 1:
				l_die_vect_cntr++;
				status_death_prob = env_1.larve_death_prob * env_1.step_size;
				break;
			case 2:
				p_die_vect_cntr++;
				status_death_prob = env_1.pupe_death_prob * env_1.step_size;
				break;
			case 3:
				if (curr_mosq->is_female == true) {
					f_die_vect_cntr++;
					status_death_prob = env_1.female_death_prob * env_1.step_size;
				}
				else {
					m_die_vect_cntr++;
					status_death_prob = env_1.male_death_prob * env_1.step_size;
				}
				break;
			}

			env_1.death_oc = gsl_ran_bernoulli(_RNG_P, status_death_prob);
			//erasing mosquitos from individual vectors if they die
			if (env_1.death_oc == 1) {
				death_cntr++;
				env_1.death_count++;
				switch (curr_mosq->status) {
				case 0:
					//cout << "\n erasing mosuito ID: " << e_die_vect_cntr;
					//cout << "\n E size"<<egg_vector.size();
					egg_vector.erase(egg_vector.begin() + e_die_vect_cntr);
					e_die_vect_cntr--;
					break;
				case 1:
					//cout << "\n erasing mosuito ID: " << l_die_vect_cntr;
					//cout << "\n L size" << larve_vector.size();
					larve_vector.erase(larve_vector.begin() + l_die_vect_cntr);
					l_die_vect_cntr--;
					break;
				case 2:
					//cout << "\n erasing mosuito ID: " << p_die_vect_cntr;
					//cout << "\n pupe size" << pupe_vector.size();
					pupe_vector.erase(pupe_vector.begin() + p_die_vect_cntr);
					p_die_vect_cntr--;
					break;
				case 3:
					if (curr_mosq->is_female == true) {
						//cout << "\n erasing mosuito ID: " << f_die_vect_cntr;
						//cout << "\n female size" << female_vector.size();
						female_vector.erase(female_vector.begin() + f_die_vect_cntr);
						f_die_vect_cntr--;
					}
					else {
						//cout << "\n erasing mosuito ID: " << m_die_vect_cntr;
						//cout << "\n male size" << male_vector.size();
						male_vector.erase(male_vector.begin() + m_die_vect_cntr);
						m_die_vect_cntr--;
						//prnt_vectr(male_vector, "Male_vector in will die");
					}
					break;
				}

				//if (curr_mosq->status == 4) {
				//	cout<<"\n erasing mosuito ID: "<<die_vect_cntr;
				//	male_vector.erase(male_vector.begin() + die_vect_cntr);
				//	die_vect_cntr--;
				//}

				//cout << "\n killing now" << curr_mosq->mosq_id << "\n";
				if (curr_mosq == head) {
					head = head->next;
					delete curr_mosq;
					curr_mosq = head;
					prev_mosq = head;
				}
				else {
					prev_mosq->next = curr_mosq->next;
					delete curr_mosq;
					curr_mosq = prev_mosq->next;
					//cout << "\n in else loop";
				}
			}
			else {
				/*if (curr_mosq->status == 4) {
					male_vector.push_back(curr_mosq);
				}*/
				prev_mosq = curr_mosq;
				curr_mosq = curr_mosq->next;
			}
			//show_mosqs();
			//cout << "\n the curr mosquito  is " << curr_mosq->mosq_id;
			//cout << "\n the temp mosquito  is " << temp->mosq_id;
			temp = curr_mosq->next; //to check if curr's next is null to execute seperate case for last mosq in list
		}//while loop

		//--------------for last mosq in list---------
		//checking death of last mosquito in the list
		//cout << "\n checking last mosq";
		//cout << "\n the mosquito in last loop before killing " << curr_mosq->mosq_id;


		//updating the indiv counters
		switch (curr_mosq->status) {
		case 0:
			e_die_vect_cntr++;
			break;
		case 1:
			l_die_vect_cntr++;
			break;
		case 2:
			p_die_vect_cntr++;
			break;
		case 3:
			if (curr_mosq->is_female == true) {
				f_die_vect_cntr++;
			}
			else {
				m_die_vect_cntr++;
			}
		}
		//die_vect_cntr++;
		total_pop = egg_vector.size() + larve_vector.size() + pupe_vector.size() + male_vector.size() + female_vector.size();
		if (total_pop <2) {
			global_zero_mosq = 1;
			cout << "\n\n\n\n zero mosqutoes alive in 1st try: \n\n\n\n";
			return;
		}

		env_1.death_oc = gsl_ran_bernoulli(_RNG_P, env_1.death_prob);
		//removing from vector for last mosq
		if (env_1.death_oc == 1) {
			switch (curr_mosq->status) {
			case 0:
				//cout << "\n erasing mosuito ID: " << e_die_vect_cntr;
				egg_vector.erase(egg_vector.begin() + e_die_vect_cntr);
				e_die_vect_cntr--;
				break;
			case 1:
				//cout << "\n erasing mosuito ID: " << l_die_vect_cntr;
				larve_vector.erase(larve_vector.begin() + l_die_vect_cntr);
				l_die_vect_cntr--;
				break;
			case 2:
				//cout << "\n erasing mosuito ID: " << p_die_vect_cntr;
				pupe_vector.erase(pupe_vector.begin() + p_die_vect_cntr);
				p_die_vect_cntr--;
				break;
			case 3:
				if (curr_mosq->is_female == true) {
					//cout << "\n erasing mosuito ID: " << f_die_vect_cntr;
					female_vector.erase(female_vector.begin() + f_die_vect_cntr);
					f_die_vect_cntr--;
				}
				else {
					//cout << "\n erasing mosuito ID: " << m_die_vect_cntr;
					male_vector.erase(male_vector.begin() + m_die_vect_cntr);
					m_die_vect_cntr--;
				}
			}
			/*if (curr_mosq->status == 4) {
				cout << "\n erasing mosuito ID: " << die_vect_cntr;
				male_vector.erase(male_vector.begin() + die_vect_cntr);
			}*/
			cout << "\n the mosquito in last loop " << curr_mosq->mosq_id;
			prev_mosq->next = NULL;
			last_mosq = prev_mosq;//adjusting last_mosq to avoid eerror while adding new mosq
			delete curr_mosq;
			cout << "\nKi;;ed last mosq in list";
			//show_mosqs();
		}

		total_pop = egg_vector.size() + larve_vector.size() + pupe_vector.size() + male_vector.size() + female_vector.size();
		if (total_pop <2) {
			global_zero_mosq = 1;
			cout << "\n\n\n\n zero mosqutoes alive in 1st try: \n\n\n\n";
			return;
		}

	}




void mosq_nxt_stg() {
	float status_next_stage_prob;
	int m_die_vect_cntr = -1;
	int f_die_vect_cntr = -1;
	int e_die_vect_cntr = -1;
	int l_die_vect_cntr = -1;
	int p_die_vect_cntr = -1;
	while (curr_mosq->next != NULL) {


		switch (curr_mosq->status) {
		case 0:
			e_die_vect_cntr++;
			status_next_stage_prob = env_1.egg_next_stg_prob* env_1.step_size;
			break;
		case 1:
			l_die_vect_cntr++;
			status_next_stage_prob = env_1.larve_next_stage_prob* env_1.step_size;
			break;
		case 2:
			p_die_vect_cntr++;
			status_next_stage_prob = env_1.pupe_next_stage_prob* env_1.step_size;
			break;
		case 3:
			if (curr_mosq->is_female == true) {
				f_die_vect_cntr++;
				status_next_stage_prob = 0;
			}
			else {
				m_die_vect_cntr++;
				status_next_stage_prob = 0;
			}
			break;
		}

		env_1.next_stg_oc = gsl_ran_bernoulli(_RNG_P, status_next_stage_prob);
		if (env_1.next_stg_oc == 1 && curr_mosq->status != 3) {
			//cout << "\n in next stage function";
				//counter for position in the vector to delete it
			switch (curr_mosq->status) {
			case 0://eggs
				curr_mosq->status = +1;
				larve_vector.push_back(curr_mosq);
				egg_vector.erase(egg_vector.begin() + e_die_vect_cntr);
				e_die_vect_cntr--;
				//checking the max number of pupa availble
				break;
			case 1://larve
				curr_mosq->status = 2;
				//cout << "\n adding pupe";
				pupe_vector.push_back(curr_mosq);
				larve_vector.erase(larve_vector.begin() + l_die_vect_cntr);
				l_die_vect_cntr--;
				break;
			case 2://pupa
				curr_mosq->status = 3;
				//cout << "\n pupa to adult from egg\n";
				if (curr_mosq->is_female == true) {
					female_vector.push_back(curr_mosq);
				}
				else {
					male_vector.push_back(curr_mosq);
				}

				pupe_vector.erase(pupe_vector.begin() + p_die_vect_cntr);
				p_die_vect_cntr--;
				break;
				//env_1.gender_oc = gsl_ran_bernoulli(_RNG_P, env_1.gender_prob);//????????????should the probability be toward male or female
				//	if (env_1.gender_oc == 1) {			//checking if the pupa will be male or female next
				//		curr_mosq->status = curr_mosq->status + 2;
				//		male_vector.push_back(curr_mosq);
				//	}
				//	else {
				//		curr_mosq->status = curr_mosq->status + 1;
				//	}
			}
		}

		curr_mosq = curr_mosq->next;
	}
	//cout << "\n the number of male mosquitos are: " << male_vector.size();
	//cout << "\n the number of female mosquitos are: " << female_vector.size();
	//cout << "\n the number of egg mosquitos are: " << egg_vector.size();
	//cout << "\n the number of lavre mosquitos are: " << larve_vector.size();
	//cout << "\n the number of pupe mosquitos are: " << pupe_vector.size();
}


void female_mating() {
	int linear_egg_flag = 0;
	double linear_egg_prob = (double)((double)1.0 - ((double)env_1.global_egg_count /(double) env_1.linear_egg_cap)); //global egg count is initally zero as there are no eggs at start of loop
	cout << "\n-------------- in female mating";
	while (curr_mosq->next != NULL) {
		if (curr_mosq->status == 3 && curr_mosq->is_female == 1 && curr_mosq->is_mated == false) { //if the mosquito is female and it is not mated
			//checking if a female mosquito is looking for a male
			env_1.find_male_oc = gsl_ran_bernoulli(_RNG_P, env_1.find_male_prob);
			if (env_1.find_male_oc == 1) {
				// if the mosquito mates sucessfully or not
				env_1.mate_oc = gsl_ran_bernoulli(_RNG_P, env_1.mate_prob* env_1.step_size);
				if (env_1.mate_oc == 1) {
					//cout << "\n prob value is : " << linear_egg_prob;
					if (env_1.linear_egg_cap > 0) { // this only works if linear egg cap is greater than zero
						if (linear_egg_prob < 0) {
							linear_egg_prob = 0;
						}
						env_1.linear_egg_oc = gsl_ran_bernoulli(_RNG_P,linear_egg_prob);
						//cout << "\n\n\n\n\n\n\n outcome is :" << env_1.linear_egg_oc;
						if (env_1.linear_egg_oc == 0) {
							//cout << "\n\n\n\n\nn\n\n\n\\n\n\n\n\\n\n inside linear egg cap";
							continue;/////////////////////////////////////work on this

						}
					}
					env_1.no_of_male = male_vector.size();
					env_1.male_rndm_mosq = gsl_rng_uniform_int(_RNG_P, env_1.no_of_male);
					if (male_vector[env_1.male_rndm_mosq]->gene_stat == 2) {
						env_1.DD_mate_cntr++;
						//treated male mosquito mating of DD(male) with RR(female)
						for (int cntr = 0; cntr < 50; cntr++) {
							add_mosq(0, 0, 1, 0, 1);
						}
						
						//cout << "\n mating with treated male mosqito\n";
						//cout << "\n female: " << curr_mosq->mosq_id << "\n mating with male: " << male_vector[env_1.male_rndm_mosq]->mosq_id << "\n";
						//unsigned long int gsl_rng_uniform_int(const gsl_rng * r, unsigned long int n)
					}

					else if (male_vector[env_1.male_rndm_mosq]->gene_stat == 1) {
						env_1.DR_mate_cntr++;
						//if the mating male is of DR(male)
						//cout << "\n inside else if of new mating\n";
						for (int cntr = 0; cntr < 25; cntr++) {
							add_mosq(0, 0, 1, 0, 1);
						}
						for (int cntr = 0; cntr < 25; cntr++) {
							add_mosq(0, 0, 1, 0, 0);
						}
					}
					else if (male_vector[env_1.male_rndm_mosq]->gene_stat == 0) {
						env_1.RR_mate_cntr++;
						//wild male mosquito mating with wild male
						//cout << "\n inside else in new mating\n";
						for (int cntr = 0; cntr < 50; cntr++) { ///changed 50 to 25 to check

							add_mosq(0, 1, 0, 0, 0); //adding 50 female
							add_mosq(0, 0, 0, 0, 0); //adding 50 male
							////env_1.gender_oc = gsl_ran_bernoulli(_RNG_P, env_1.gender_prob); // gender to the egg
							//if (env_1.gender_oc == 1) {//if gender is female
							//	add_mosq(0, 1, 0, 0,0);
							//}
							//else {
							//	add_mosq(0, 0, 0, 0, 0);
							//}


						}
						//cout << "\n mated with wild male mosqito8888888888888888888888888888888888888888888888888888888\n";
						//cout << "\n female: " << curr_mosq->mosq_id << "\n mating with male: " << male_vector[env_1.male_rndm_mosq]->mosq_id << "\n";
						//show_mosqs();
					}
					curr_mosq->is_mated = 1;
					male_vector[env_1.male_rndm_mosq]->is_mated = 1;
				}
			}
		}
		curr_mosq = curr_mosq->next;
	}

}

void larve_id_coll_and_shuff() {
	curr_mosq = head;
	prev_mosq = head;
	larve_ids.clear();
	while (curr_mosq->next != NULL) {
		if (curr_mosq->status == 1) {
			larve_ids.push_back(curr_mosq->mosq_id);
		}
		prev_mosq = curr_mosq;
		curr_mosq = curr_mosq->next;
	}

	//collected larve_ids from all mosq list.

	//creating array for corresponding larve vector

	//int arr_size = larve_ids.size();
	//int* larve_id_array = larve_ids.data();
	//
	//for (int i = 0; i < arr_size; i++) {
	//	cout << " " << larve_id_array[i];
	//}

	//gsl_ran_shuffle(_RNG_P,  larve_id_array, arr_size, sizeof(int));

	//cout << "\n\n after shuffle\n";
	//for (int i = 0; i < arr_size; i++) {
	//	cout << " " << larve_id_array[i];
	//}
}

void egg_id_coll_and_shuff() {
	curr_mosq = head;
	prev_mosq = head;
	egg_ids.clear();
	eggs_female_vector.clear();
	eggs_male_DR_vector.clear();
	eggs_male_RR_vector.clear();

	while (curr_mosq->next != NULL) {
		if (curr_mosq->status == 0) {
			egg_ids.push_back(curr_mosq->mosq_id);

			//code for collection gene and gender of eggs
			if (curr_mosq->is_female == 1 && curr_mosq->gene_stat == 0) {//for female eggs
				eggs_female_vector.push_back(curr_mosq->mosq_id);

			}

			//code for coll of male RR eggs
			else if (curr_mosq->is_female == 0 && curr_mosq->gene_stat == 0) {
				eggs_male_RR_vector.push_back(curr_mosq->mosq_id);
			}

			//code for coll of male DR eggs
			else if (curr_mosq->is_female == 0 && curr_mosq->gene_stat == 1) {
				eggs_male_DR_vector.push_back(curr_mosq->mosq_id);
			}
			
		}
		prev_mosq = curr_mosq;
		curr_mosq = curr_mosq->next;
	}
}

void kill_excess_larve() {
	int l_die_vect_cntr = -1;
	curr_mosq = head;
	prev_mosq = head;
	int larve_counter = -1; // to keep track for iteration in loop to start killing
	int total_pop = egg_vector.size() + larve_vector.size() + pupe_vector.size() + male_vector.size() + female_vector.size(); //finding the population of all mosquitoes
	cout << "\n the total population of all mosquitoes is: " << total_pop;
	//for (int i = 0; i < total_pop; i++)
	//prnt_vectr(larve_vector, "Larve");


	int arr_size = larve_ids.size();
	int* larve_id_array = larve_ids.data();
	int mosq_id_check = 0;
	//cout << "\n before shuffle";
	//for (int i = 0; i < arr_size; i++) {//after shuffle of larve
	//	cout << " " << larve_id_array[i];
	//}

	gsl_ran_shuffle(_RNG_P, larve_id_array, arr_size, sizeof(int));
	//cout << "\n after suffle\n";
//	for (int i = 0; i < arr_size; i++) {//after shuffle of larve
//
//		if (i > 100) {
//			cout << "--";
//		}
//	cout << " " << larve_id_array[i];
//}

	while (curr_mosq->next != NULL) {

		if (curr_mosq->status == 1) {

			for (int i = 0; i < env_1.max_larve_cap; i++) {
				if (curr_mosq->mosq_id == larve_id_array[i]) { //safe mosq
					//cout << "\n " << curr_mosq->mosq_id;
					mosq_id_check = 1;
					break;
				}
			}

			l_die_vect_cntr++;

			if (mosq_id_check == 1) {//mosq is in safe list
				mosq_id_check = 0;
				//cout << "\n not killing " << curr_mosq->mosq_id;
				prev_mosq = curr_mosq;
				curr_mosq = curr_mosq->next;
			}
			else {
				//cout << "\n killing not safe mosq" << curr_mosq->mosq_id;
				if (curr_mosq == head) {
					head = head->next;
					delete curr_mosq;
					curr_mosq = head;
					prev_mosq = head;
				}
				else {
					prev_mosq->next = curr_mosq->next;
					delete curr_mosq;
					curr_mosq = prev_mosq->next;
					//cout << "\n in else loop of kill _L";
				}

				larve_vector.erase(larve_vector.begin() + l_die_vect_cntr);
				//cout << "\n current vector pos: " << l_die_vect_cntr;
				l_die_vect_cntr--;
			}

			}

		else {
			prev_mosq = curr_mosq;
			curr_mosq = curr_mosq->next;
		}

	}
}


void kill_excess_eggs() {
	int e_die_vect_cntr = -1;
	curr_mosq = head;
	prev_mosq = head;
	int egg_counter = -1; // to keep track for iteration in loop to start killing
	//int total_pop = egg_vector.size() + larve_vector.size() + pupe_vector.size() + male_vector.size() + female_vector.size(); //finding the population of all mosquitoes
	//cout << "\n the total population of all mosquitoes is: " << total_pop;
	////for (int i = 0; i < total_pop; i++)
	//prnt_vectr(egg_vector, "egg");

	int female_egg_arr_size = eggs_female_vector.size();
	int male_RR_eggs_arr_size = eggs_male_RR_vector.size();
	int male_DR_eggs_arr_size = eggs_male_DR_vector.size();

	int* female_egg_id_array = eggs_female_vector.data();
	int* male_RR_egg_id_array = eggs_male_RR_vector.data();
	int* male_DR_egg_id_array = eggs_male_DR_vector.data();

	int arr_size = egg_ids.size();
	int* egg_id_array = egg_ids.data();
	int mosq_id_check = 0;
	//cout << "\nbefore shuffle\n";

	//cout << "\n\n\n below are female eggs\n";

	//for (int i = 0; i < female_egg_arr_size; i++) {//after shuffle of larve
	//	cout << " " << female_egg_id_array[i];
	//}

	//cout << "\n\n\n below are RR male egs\n";

	//for (int i = 0; i < male_RR_eggs_arr_size; i++) {//after shuffle of larve
	//	cout << " " << male_RR_egg_id_array[i];
	//}

	//cout << "\n\n\n below are DR male egs\n";
	//for (int i = 0; i < male_DR_eggs_arr_size; i++) {//after shuffle of larve
	//	cout << " " << male_DR_egg_id_array[i];
	//}

	//for (int i = 0; i < arr_size; i++) {//after shuffle of larve
	//	cout << " " << egg_id_array[i];
	//}

	cout << "\n the proportions are Female: "<<eggs_female_vector.size()<<"::"<<female_egg_arr_size << " DR: " << eggs_male_DR_vector.size()<<"::"<<male_DR_eggs_arr_size << " RR: " << eggs_male_RR_vector.size() <<"::"<<male_RR_eggs_arr_size << " total mosqs: " << egg_ids.size();
	cout << "\n\n before shuffle\n";
	cout << "\n the proprtions are: female " << (double)female_egg_arr_size / (double)arr_size;
	cout<<"\nmale DR "<<(double)male_DR_eggs_arr_size/ (double)arr_size;
	cout<<"\n Male RR "<< (double)male_RR_eggs_arr_size / (double)arr_size;

	gsl_ran_shuffle(_RNG_P, egg_id_array, arr_size, sizeof(int));

	//for (int i = 0; i < arr_size; i++) {//after shuffle of larve
	//	cout << " " << egg_id_array[i];
	//}
	while (curr_mosq->next != NULL) {

		if (curr_mosq->status == 0) {

			for (int i = 0; i < env_1.max_egg_cap; i++) {
				if (curr_mosq->mosq_id == egg_id_array[i]) { //safe mosq
					mosq_id_check = 1;
					break;
				}
			}

			e_die_vect_cntr++;

			if (mosq_id_check == 1) {//mosq is in safe list
				//cout << "\n not killing egg" << curr_mosq->mosq_id;
				mosq_id_check = 0;
				prev_mosq = curr_mosq;
				curr_mosq = curr_mosq->next;
			}
			else {
				//cout << "\n killing egg " << curr_mosq->mosq_id;
				if (curr_mosq == head) {
					head = head->next;
					delete curr_mosq;
					curr_mosq = head;
					prev_mosq = head;
				}
				else {
					prev_mosq->next = curr_mosq->next;
					delete curr_mosq;
					curr_mosq = prev_mosq->next;
					//cout << "\n in else loop of kill _L";
				}

				egg_vector.erase(egg_vector.begin() + e_die_vect_cntr);
				//cout << "\n current vector pos: " << l_die_vect_cntr;
				e_die_vect_cntr--;
			}

		}

		else {
			prev_mosq = curr_mosq;
			curr_mosq = curr_mosq->next;
		}

	}
	egg_id_coll_and_shuff();

	int female_egg_arr_size1 = eggs_female_vector.size();
	int male_RR_eggs_arr_size1 = eggs_male_RR_vector.size();
	int male_DR_eggs_arr_size1 = eggs_male_DR_vector.size();

	int* female_egg_id_array1 = eggs_female_vector.data();
	int* male_RR_egg_id_array1= eggs_male_RR_vector.data();
	int* male_DR_egg_id_array1 = eggs_male_DR_vector.data();

	int arr_size1 = egg_ids.size();
	int* egg_id_array1 = egg_ids.data();

	//cout << "\nafter shuffle\n";
	//cout << "\n below are female eggs\n";
	//for (int i = 0; i < female_egg_arr_size1; i++) {//after shuffle of larve
	//	cout << " " << female_egg_id_array1[i];
	//}

	//cout << "\n\n\n below are RR male eggs\n";

	//for (int i = 0; i < male_RR_eggs_arr_size1; i++) {//after shuffle of larve
	//	cout << " " << male_RR_egg_id_array1[i];
	//}

	//cout << "\n\n\n below are DR male egs\n";
	//for (int i = 0; i < male_DR_eggs_arr_size1; i++) {//after shuffle of larve
	//	cout << " " << male_DR_egg_id_array1[i];
	//}

	//for (int i = 0; i < arr_size; i++) {//after shuffle of larve
	//	cout << " " << egg_id_array[i];
	//}

	cout << "\n the proportions are Female: " << eggs_female_vector.size() << "::" << female_egg_arr_size1 << " DR: " << eggs_male_DR_vector.size() << "::" << male_DR_eggs_arr_size1 << " RR: " << eggs_male_RR_vector.size() << "::" << male_RR_eggs_arr_size1 << " total mosqs: " << egg_ids.size();
	cout << "\n after shuffle\n";
	cout << "\n the proprtions are: female " << (double)female_egg_arr_size1 / (double)arr_size1;
	cout << "\nmale DR " << (double)male_DR_eggs_arr_size1 / (double)arr_size1;
	cout << "\n Male RR " << (double)male_RR_eggs_arr_size1 / (double)arr_size1;




}
		


void kill_excess_rndm_eggs() {// function to kill random eggs based on proportion of male RR and RD amd female

	//adult=3,0=male,is_treated=0, is_mated=0, gene_stat=0 RR /1 DR




	int e_die_vect_cntr = -1;
	curr_mosq = head;
	prev_mosq = head;
	int egg_counter = -1; // to keep track for iteration in loop to start killing
	//int total_pop = egg_vector.size() + larve_vector.size() + pupe_vector.size() + male_vector.size() + female_vector.size(); //finding the population of all mosquitoes
	//cout << "\n the total population of all mosquitoes is: " << total_pop;
	////for (int i = 0; i < total_pop; i++)
	//prnt_vectr(egg_vector, "egg");


	int arr_size = egg_ids.size();
	int* egg_id_array = egg_ids.data();
	int mosq_id_check = 0;
	cout << "\nbefore shuffle\n";
	//for (int i = 0; i < arr_size; i++) {//after shuffle of larve
	//	cout << " " << egg_id_array[i];
	//}

	gsl_ran_shuffle(_RNG_P, egg_id_array, arr_size, sizeof(int));

	cout << "\nafter shuffle\n";
	//for (int i = 0; i < arr_size; i++) {//after shuffle of larve
	//	cout << " " << egg_id_array[i];
	//}
	while (curr_mosq->next != NULL) {

		if (curr_mosq->status == 0) {

			for (int i = 0; i < env_1.max_egg_cap; i++) {
				if (curr_mosq->mosq_id == egg_id_array[i]) { //safe mosq
					mosq_id_check = 1;
					break;
				}
			}

			e_die_vect_cntr++;

			if (mosq_id_check == 1) {//mosq is in safe list
				//cout << "\n not killing egg" << curr_mosq->mosq_id;
				mosq_id_check = 0;
				prev_mosq = curr_mosq;
				curr_mosq = curr_mosq->next;
			}
			else {
				//cout << "\n killing egg " << curr_mosq->mosq_id;
				if (curr_mosq == head) {
					head = head->next;
					delete curr_mosq;
					curr_mosq = head;
					prev_mosq = head;
				}
				else {
					prev_mosq->next = curr_mosq->next;
					delete curr_mosq;
					curr_mosq = prev_mosq->next;
					//cout << "\n in else loop of kill _L";
				}

				egg_vector.erase(egg_vector.begin() + e_die_vect_cntr);
				//cout << "\n current vector pos: " << l_die_vect_cntr;
				e_die_vect_cntr--;
			}

		}

		else {
			prev_mosq = curr_mosq;
			curr_mosq = curr_mosq->next;
		}

	}
}








void add_treated_male(int mosq_count) {
	for (int i = 0; i < mosq_count; i++) {
		add_mosq(3, 0, 1, 0, 2);
	}
}

bool comparePtrToNode(mosq* a, mosq* b) { return (a->mosq_id < b->mosq_id); }



/////////////
int main() {

	string file_names_list = "80_percent_sim_";
	string temp_file_names;
	

	//for (int sim_cnt; sim_cnt < 50; sim_cnt++) {



		alloc_and_setup_gsl_random_generator();
		seed_gsl_random_generator(0);

		cout << "the intialized count :";
		curr_mosq = head;
		prev_mosq = head;
		int day_cntr = 0;
		//int Male_DD = 0, Male_DR = 0, Male_RR = 0;

		double divs;
		read_variables();
		cout << "\n the initial mosqs count is :" << env_1.mate_prob;
		cout << "\n" << env_1.initial_mosqs;
		cout << "\n" << env_1.max_larve_cap;
		cout << "\n" << env_1.mate_prob;
		cout << "\n" << env_1.treated_prob;
		cout << "\n" << env_1.female_death_prob;
		cout << "\n" << env_1.male_death_prob;
		cout << "\n" << env_1.egg_death_prob;
		cout << "\n" << env_1.larve_death_prob;
		cout << "\n" << env_1.pupe_death_prob;
		cout << "\n" << env_1.gender_prob;
		cout << "\n" << env_1.find_male_prob;
		cout << "\n" << env_1.egg_next_stg_prob;
		cout << "\n" << env_1.larve_next_stage_prob;
		cout << "\n" << env_1.pupe_next_stage_prob;


		for (int i = 0; i < env_1.initial_mosqs; i++) {
			//k = gsl_ran_bernoulli(_RNG_P, 0.50);  // stoping this as we are bringing all the mosquitos to a new environment
			if (env_1.k == 1) {
				env_1.gender_oc = gsl_ran_bernoulli(_RNG_P, env_1.gender_prob);
				env_1.mate_oc = gsl_ran_bernoulli(_RNG_P, env_1.initial_mate_prob);
				env_1.treated_oc = gsl_ran_bernoulli(_RNG_P, env_1.treated_prob);
				if (env_1.gender_oc == 1) {
					add_mosq(3, 0, 0, env_1.mate_oc, 0);
				}
				else {
					add_mosq(3, 1, 0, env_1.mate_oc, 0);
				}
			}
		}

		cout << "\n Male vector size" << male_vector.size();

		//add_mosq(3, 0, 0, env_1.mate_oc,/* 1);
		//add_mosq(3, 0, 0, env_1.mate_oc, 2);*/
		//add_mosq(3, 0, 0, env_1.mate_oc, 0);



		//to calculate seperate count for DD, DR,RR
		cout << "\n Male vector size" << male_vector.size();
		Male_DD.clear();
		Male_DR.clear();
		Male_RR.clear();
		for (int si = 0; si < male_vector.size(); si++) {
			cout << "\n" << male_vector[si]->mosq_id << "  ----  " << male_vector[si]->gene_stat;
			if (male_vector[si]->gene_stat == 2) {
				//Male_DD++;
				Male_DD.push_back(male_vector[si]->mosq_id);
			}
			else if (male_vector[si]->gene_stat == 1) {
				Male_DR.push_back(male_vector[si]->mosq_id);
			}
			else if (male_vector[si]->gene_stat == 0) {
				Male_RR.push_back(male_vector[si]->mosq_id);
			}

		}

		Male_count_DD.push_back(Male_DD.size());
		Male_count_DR.push_back(Male_DR.size());
		Male_count_RR.push_back(Male_RR.size());

		//Male_DD.clear();
		//Male_DR.clear();
		//Male_RR.clear();
		//for (int si = 0; si < male_vector.size(); si++) {
		//	cout <<"\n" << male_vector[si]->mosq_id <<"  ----  "<< male_vector[si]->gene_stat;
		//	if (male_vector[si]->gene_stat == 1) {
		//		//Male_DD++;
		//		
		//		Male_DD.push_back(male_vector[si]->mosq_id);
		//	}
		//	else if (male_vector[si]->gene_stat == 2) {
		//		Male_DD.push_back(male_vector[si]->mosq_id);
		//	}
		//	else if (male_vector[si]->gene_stat == 0) {
		//		Male_DD.push_back(male_vector[si]->mosq_id);
		//	}

		//}

		int males_add_after_equi;

		cout << "\n count is";
		cout << Male_DD.size() << "\n" << Male_DR.size() << "\n" << Male_RR.size();


		/*show_mosqs();
		prnt_all_vect();*/

		/*cout << "\n male mosquitoes are:";
		for (auto it = male_vector.begin(); it != male_vector.end(); it++) {
			cout << (*it)->mosq_id << "\n";
		}*/
		cout << "\n===============" << mosq_count;
		//loooping over time
		for (int day = 0; day < env_1.num_of_days_to_expmnt; day++) {

			for (int time = 0; time < env_1.num_of_times_in_each_day; time++) {



				cout << "\n loop counter " << day << " ------- " << time;
				cout << "\n male vector size : " << male_vector.size();
				int i = 0;
				temp = head;
				curr_mosq = head;
				prev_mosq = head;
				cout << "\n-----------killing started------------\n";
				cout << "\n male vector size : " << male_vector.size();
				will_mosq_die();//checking each mosquito survives each time loop and killing them based on probability

				if (global_zero_mosq == 1) {// try catch throw : when there are zero mosq, this will exit the loop and send data to files.
					cout << "\n exiting the inner ---for--- loop to write data to files\n";
					break;
				}

				prnt_all_vect();
				//adding l* implementation for larve to kill excess



				temp = head;
				i = 0;
				temp = head;
				curr_mosq = head;
				prev_mosq = head;
				//prnt_all_vect();
				mosq_nxt_stg();// checkimg if mosquito forwards to next stage like from egg to larve
				std::sort(male_vector.begin(), male_vector.end(), comparePtrToNode);
				std::sort(female_vector.begin(), female_vector.end(), comparePtrToNode);
				std::sort(egg_vector.begin(), egg_vector.end(), comparePtrToNode);
				std::sort(larve_vector.begin(), larve_vector.end(), comparePtrToNode);
				std::sort(pupe_vector.begin(), pupe_vector.end(), comparePtrToNode);
				prnt_all_vect();

				// for all mosquitos
				curr_mosq = head;

				if (male_vector.size() > 0 && female_vector.size() > 0) {
					female_mating();// checking if female mates with male and lays eggs
				}
				else {
					cout << "\n\n\n\n\n\n no adult mosqitoes left for mating\n";
				}

				prnt_all_vect();

				if (egg_vector.size() > env_1.max_egg_cap) {
					egg_id_coll_and_shuff();
					kill_excess_eggs();
				}


				cout << "\n curr larve count: " << larve_vector.size();
				if (larve_vector.size() > env_1.max_larve_cap) {
					larve_id_coll_and_shuff();
					kill_excess_larve();
				}

				cout << "\n the mosquitos after deletion --------------------------- far are \n";
				/*		cout << "\n Mosquito ID: ";
						cout << " Mosquito stage: ";
						cout << " treated status: ";
						cout << " mated status: ";*/
				show_mosqs();

				//cout << "\n\n the number of male mosquitos are :" << male_vector.size() << endl;
			/*	cout << "\n" << "the male mosquitos left are \n";
				for (auto it = male_vector.begin(); it != male_vector.end(); it++) {
					cout << *it << "\t" << (*it)->mosq_id << "\n";
				}*/

				//cout << "\n the number of male mosquitos are: " << male_vector.size();
				//cout << "\n the number of female mosquitos are: " << female_vector.size();
				//cout << "\n the number of egg mosquitos are: " << egg_vector.size();
				//cout << "\n the number of lavre mosquitos are: " << larve_vector.size();
				//cout << "\n the number of pupe mosquitos are: " << pupe_vector.size();


				/*if (female_vector.size() < 10) {
					cout << "\n female vector size is low";
					cout << female_vector.size();
					cin >> pause_for_interupt;
				}*/

				env_1.global_egg_count = egg_vector.size(); // to keep track of egg count after each time step.

			}





			if (global_zero_mosq == 1) {// try catch throw : when there are zero mosq, this will exit the loop and send data to files.
				cout << "\n exiting the outer ---for--- loop to write data to files\n";
				break;
			}

			if (day == env_1.equlibrium_days) {
				env_1.adult_males_on_equibrilium = male_vector.size();
				males_add_after_equi = (int)((double)env_1.adult_males_on_equibrilium * (double)env_1.percent_males_DD_add * (double)0.01);

			}

			//to add mosquitoes in particualr interval

			if (day % env_1.intrvl_of_days_to_add_treated_mosq == 0 && day > env_1.equlibrium_days) {
				cout << "\n male vector size : " << male_vector.size();
				cout << "\n ------------adding treated mosqs";
				//add_treated_male(env_1.num_of_mosq_to_add);

				for (int i = 0; i < males_add_after_equi; i++) {
					add_mosq(3, 0, 0, 0, 2);
				}

			}
			cout << "\n male vector size : " << male_vector.size();

			male_count.push_back(male_vector.size());
			female_count.push_back(female_vector.size());
			egg_count.push_back(egg_vector.size());
			larve_count.push_back(larve_vector.size());
			pupe_count.push_back(pupe_vector.size());
			time_count.push_back(day_cntr++);



			female_ratio.push_back((double)female_vector.size() / (double)larve_vector.size());
			male_ratio.push_back((double)male_vector.size() / (double)larve_vector.size());
			larve_ratio.push_back((double)larve_vector.size() / (double)larve_vector.size());
			pupe_ratio.push_back((double)pupe_vector.size() / (double)larve_vector.size());

			Male_DD.clear();
			Male_DR.clear();
			Male_RR.clear();
			for (int si = 0; si < male_vector.size(); si++) {
				if (male_vector[si]->gene_stat == 2) {
					Male_DD.push_back(male_vector[si]->mosq_id);
				}
				else if (male_vector[si]->gene_stat == 1) {
					Male_DR.push_back(male_vector[si]->mosq_id);
				}
				else if (male_vector[si]->gene_stat == 0) {
					Male_RR.push_back(male_vector[si]->mosq_id);
				}

			}

			Male_count_DD.push_back(Male_DD.size());
			Male_count_DR.push_back(Male_DR.size());
			Male_count_RR.push_back(Male_RR.size());

			cout << "\n count is";
			cout << Male_DD.size() << "\n" << Male_DR.size() << "\n" << Male_RR.size();


			//mate counting for each day for DD,DR,RR with female
			DD_mate_count_vect.push_back(env_1.DD_mate_cntr);
			env_1.DD_mate_cntr = 0;
			DR_mate_count_vect.push_back(env_1.DR_mate_cntr);
			env_1.DR_mate_cntr = 0;
			RR_mate_count_vect.push_back(env_1.RR_mate_cntr);
			env_1.RR_mate_cntr = 0;


		}

		// add function to seperate the read and write to files
		//Male_DD.clear();
		//Male_DR.clear();
		//Male_RR.clear();
		//for (int si = 0; si < male_vector.size(); si++) {
		//	if (male_vector[si]->gene_stat == 2) {
		//		//Male_DD++;
		//		Male_DD.push_back(male_vector[si]->mosq_id);
		//	}
		//	else if (male_vector[si]->gene_stat == 1) {
		//		Male_DD.push_back(male_vector[si]->mosq_id);
		//	}
		//	else if (male_vector[si]->gene_stat == 0) {
		//		Male_DD.push_back(male_vector[si]->mosq_id);
		//	}

		//}

		//int m_DD = Male_DD;
		//int m_DR=
		cout << "\n files in write";
		std::vector<std::pair<std::string, std::vector<int>>> vals3 = { {"time",time_count }, {"Egg", egg_count}, {"Larve", larve_count},{"Pupe",pupe_count},{"Adult_male",male_count},{"Adult_female",female_count},{"MALE_DD",Male_count_DD},{"MALE_DR",Male_count_DR},{"MALE_RR",Male_count_RR},{"DD_mate_count",DD_mate_count_vect},{"DR_mate_count",DR_mate_count_vect},{"RR_mate_count",RR_mate_count_vect} };
		std::vector<std::pair<std::string, std::vector<double>>> ratios = { { "Larve_ratio",larve_ratio }, { "Pupe_ratio",pupe_ratio}, { "Male_ratio",male_ratio}, { "Female_ratio",female_ratio} };

		cout << "\n number of adults on equibrilium day: " << env_1.adult_males_on_equibrilium << "\n added with percent " << males_add_after_equi;
		cout << "\n files in write";
		temp_file_names = file_names_list + "_31" + ".csv";
		write_csv(temp_file_names, vals3);
		file_counter++;

		cout << "\n files in write";
		//write_double_CSV(larve_ratio, "data.csv");
		//write_csv("ratios.csv", ratios);
		write_d_csv("data1.csv", ratios);
		cout << "\n files in writing completed";
		male_count.clear();
		male_count.erase(male_count.begin(), male_count.end());
		male_count.shrink_to_fit();

		female_count.clear();
		female_count.erase(female_count.begin(), female_count.end());
		female_count.shrink_to_fit();

		egg_count.clear();
		egg_count.erase(egg_count.begin(), egg_count.end());
		egg_count.shrink_to_fit();

		larve_count.clear();
		larve_count.erase(larve_count.begin(), larve_count.end());
		larve_count.shrink_to_fit();

		pupe_count.clear();
		pupe_count.erase(pupe_count.begin(), pupe_count.end());
		pupe_count.shrink_to_fit();

		time_count.clear();
		time_count.erase(time_count.begin(), time_count.end());
		time_count.shrink_to_fit();


		male_vector.clear();
		male_vector.erase(male_vector.begin(), male_vector.end());
		male_vector.shrink_to_fit();

		larve_vector.clear();
		larve_vector.erase(larve_vector.begin(), larve_vector.end());
		larve_vector.shrink_to_fit();

		egg_vector.clear();
		egg_vector.erase(egg_vector.begin(), egg_vector.end());
		egg_vector.shrink_to_fit();

		pupe_vector.clear();
		pupe_vector.erase(pupe_vector.begin(), pupe_vector.end());
		pupe_vector.shrink_to_fit();

		female_vector.clear();
		female_vector.erase(female_vector.begin(), female_vector.end());
		female_vector.shrink_to_fit();

		female_ratio.clear();
		female_ratio.erase(female_ratio.begin(), female_ratio.end());
		female_ratio.shrink_to_fit();

		male_ratio.clear();
		male_ratio.erase(male_ratio.begin(), male_ratio.end());
		male_ratio.shrink_to_fit();

		larve_ratio.clear();
		larve_ratio.erase(larve_ratio.begin(), larve_ratio.end());
		larve_ratio.shrink_to_fit();

		pupe_ratio.clear();
		pupe_ratio.erase(pupe_ratio.begin(), pupe_ratio.end());
		pupe_ratio.shrink_to_fit();


		Male_DD.clear();
		Male_DD.erase(Male_DD.begin(), Male_DD.end());
		Male_DD.shrink_to_fit();

		Male_DR.clear();
		Male_DR.erase(Male_DR.begin(), Male_DR.end());
		Male_DR.shrink_to_fit();

		Male_RR.clear();
		Male_RR.erase(Male_RR.begin(), Male_RR.end());
		Male_RR.shrink_to_fit();

		Male_count_DD.clear();
		Male_count_DD.erase(Male_count_DD.begin(), Male_count_DD.end());
		Male_count_DD.shrink_to_fit();

		Male_count_DR.clear();
		Male_count_DR.erase(Male_count_DR.begin(), Male_count_DR.end());
		Male_count_DR.shrink_to_fit();

		Male_count_RR.clear();
		Male_count_RR.erase(Male_count_RR.begin(), Male_count_RR.end());
		Male_count_RR.shrink_to_fit();

		larve_ids.clear();
		larve_ids.erase(larve_ids.begin(), larve_ids.end());
		larve_ids.shrink_to_fit();


		egg_ids.clear();
		egg_ids.erase(egg_ids.begin(), egg_ids.end());
		egg_ids.shrink_to_fit();

		eggs_female_vector.clear();
		eggs_female_vector.erase(eggs_female_vector.begin(), eggs_female_vector.end());
		eggs_female_vector.shrink_to_fit();

		eggs_male_RR_vector.clear();
		eggs_male_RR_vector.erase(eggs_male_RR_vector.begin(), eggs_male_RR_vector.end());
		eggs_male_RR_vector.shrink_to_fit();

		eggs_male_DR_vector.clear();
		eggs_male_DR_vector.erase(eggs_male_DR_vector.begin(), eggs_male_DR_vector.end());
		eggs_male_DR_vector.shrink_to_fit();

		DD_mate_count_vect.clear();
		DD_mate_count_vect.erase(DD_mate_count_vect.begin(), DD_mate_count_vect.end());
		DD_mate_count_vect.shrink_to_fit();

		DR_mate_count_vect.clear();
		DR_mate_count_vect.erase(DR_mate_count_vect.begin(), DR_mate_count_vect.end());
		DR_mate_count_vect.shrink_to_fit();

		RR_mate_count_vect.clear();
		RR_mate_count_vect.erase(RR_mate_count_vect.begin(), RR_mate_count_vect.end());
		RR_mate_count_vect.shrink_to_fit();



		curr_mosq = head;

		cout << "\n " << head->mosq_id;
		cout << "\n " << curr_mosq->mosq_id;
		cout << "\n " << last_mosq->mosq_id;
	//}

}

