extern "C" {
	#include <unistd.h>
	#include <sys/types.h>
	#include <sys/wait.h>
	#include <pthread.h>
}

#include <climits>
#include <iostream>
#include <map>
#include <unordered_map>
#include <cstdint>
#include <unordered_set>
#include <fstream>
#include <ctime>
#include <random>
#include <cmath>
#include <queue>

#include "../fast-cpp-csv-parser/csv.h"

#include "settings.h"

using namespace std;

//nasty global variable, do not write after init
struct timespec start_time;
const int clock_status = clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);


class timetable_t {
private:
	uint16_t ***data;
	int space = 300;
	static const uint16_t internal_empty_value = UINT16_MAX;//should be the biggest value - most expensive flight

	int get_new_space() const
	{
		//TODO: align size with expected problem size
		return space*2;
	}

	timetable_t (const timetable_t &obj) {//private copy constructor, just in case
		data = NULL;
		throw obj;
	}
	
public:
	inline static unsigned get_empty_value() {
		return UINT_MAX;
	}

	timetable_t() {
		int i, j, k;
		//no checking, let it crash ;-)
		data = (uint16_t***)malloc(space*sizeof(**data));
		for (i=0; i<space; ++i) {
			data[i] = (uint16_t**)malloc(space*sizeof(*data));
			for (j=0; j<space; ++j) {
				data[i][j] = (uint16_t*)malloc(space*sizeof(data));
				for (k=0; k<space; ++k) {
					data[i][j][k]=internal_empty_value;
				}
			}
		}
	}

	~timetable_t() {
		int i, j;
		for (i=0; i<space; ++i) {
			for (j=0; j<space; ++j) {
				free(data[i][j]);
			}
			free(data[i]);
		}
		free(data);
	}


	inline unsigned get_price(int from, int to, int day) const {
//#ifdef DEBUG
//		return data.at(from).at(day).at(to);
//#else
		const uint16_t raw_data = data[from][day][to];
		if (raw_data == internal_empty_value) {
			return get_empty_value();
		}
		return raw_data;//type-casted to int
//#endif
	}


	void set_price(int from, int to, int day, uint16_t price) {
		//TODO:resizing segfaults
		if ((from >= space) || (to >= space)) {
			int i, j, k;
			int old_space = space;
			space = get_new_space();
			data = (uint16_t***)realloc(data, space*sizeof(**data));
			for (i=old_space; i<space; ++i) {
				data[i] = (uint16_t**)malloc(space*sizeof(*data));
				for (j=0; j<space; ++j) {
					data[i][j] = (uint16_t*)malloc(space*sizeof(data));
					for (k=0; k<space; ++k) {
						data[i][j][k]=internal_empty_value;
					}
				}
			}
			for (i=0; i<space; ++i) {
				data[i] = (uint16_t**)realloc(data[i], space*sizeof(*data));
				for (j=old_space; j<space; ++j) {
					data[i][j] = (uint16_t*)malloc(space*sizeof(data));
					for (k=0; k<space; ++k) {
						data[i][j][k]=internal_empty_value;
					}
				}
				for (j=0; j<space; ++j) {
					data[i][j] = (uint16_t*)realloc(data[i][j], space*sizeof(data));
					for (k=old_space; k<space; ++k) {
						data[i][j][k]=internal_empty_value;
					}
				}
			}
		}
		data[from][day][to] = price;
	}
};


class problem_t {
public:
	string starting_city_string;
	int starting_city_index;
	
	string get_city_name(int index) const {
		for (auto iter=city_names.begin(); iter != city_names.end(); ++iter) {
			if (iter->second == index) {
				return iter->first;
			}
		}
		cerr << "ERROR: city not found";
		exit(-1);
	}


	problem_t() {
		cin >> starting_city_string;
		//city_names[starting_city_string] = 0;
		while (cin.get()!='\n') {//do not break csv by starting with a stray newline
#ifdef DEBUG
			if (!cin) {
				cerr << "ERROR: could not read stdin\n";
				exit(-1);
			}
#endif
		}
		//io::LineReader in("stdin", cin);
		//while(char*line = in.next_line()){
		//	cout << line;
		//}		
#ifdef DEBUG
	#define CSV_OVERFLOW throw_on_overflow
#else
	#define CSV_OVERFLOW ignore_overflow
#endif
		io::CSVReader<4, io::trim_chars<' '>, io::no_quote_escape<' '>, io::CSV_OVERFLOW, io::no_comment> in("stdin", cin);
		in.set_header("from", "to", "day", "price");
		string from, to; unsigned short day=0, price=0;
		//int sum=0;
		while(in.read_row(from, to, day, price)){
			//cout << from << " " << to << " " << day << " " << price << endl;
			//if (from[0] || to[0] || day || price) {
			//	sum++;
			//}
			
			//if (this->city_names.count(from) == 0 ) {
			//	this->city_names[from] = this->city_names.size();
			//}
			
			//if (this->city_names.count(to) == 0 ) {
			//cout << "pridavam "<<city_names.size() <<endl;
			//city_names[to] = city_names.size();
			//city_names.emplace(from, city_names.size());
			//city_names.emplace(to, city_names.size());
			city_names.insert(pair<string, int>(from, city_names.size()));
			city_names.insert(pair<string, int>(to, city_names.size()));

			timetable.set_price(city_names[from], city_names[to], day, price);
			//cout << "pridane "<<city_names.size() <<endl;
			//}
		}
		//for (auto iter=city_names.begin(); iter != city_names.end(); ++iter) {
		//	cout << iter->first << " " << iter->second << endl;
		//}

		//cout<<sum;
		starting_city_index = city_names[starting_city_string];
	}

	int problem_size() const {
		return city_names.size();
	}

	/*void print_table() {
		timetable.print_table(city_names);
	}*/

	inline unsigned get_price(int from, int to, int day) const {
		return timetable.get_price(from, to, day);
	}
	
	inline static unsigned get_empty_value() {
		return timetable_t::get_empty_value();
	}

		
private:
	timetable_t timetable;
	unordered_map<string, int> city_names;
	//int city_count = 0;
	
	problem_t(const problem_t &obj) {//private copy constructor
		throw obj;
	}

};

inline void un_set_remove(unordered_set<int> &set, const int elem)
{
	//cout << "pred mazanim: "<<set.size();
	set.erase(set.find(elem));
	//cout << " po mazani: "<<set.size()<<endl;
}

unsigned get_path_cost(const vector<int> &path, const problem_t &problem) {
	auto iter = path.begin();
	int prev_city = *iter;
	int day = 0;
	unsigned sum = 0;
	for (++iter; iter != path.end(); ++iter) {
		unsigned cost = problem.get_price(prev_city, *iter, day);
		if (cost == problem.get_empty_value()) {
			return  problem.get_empty_value();
		}
		sum += cost;
		++day;
		prev_city = *iter;
	}
	return sum;
}


class solution_t {
	public:
	vector<int> path;
	bool valid = false;
	unsigned cost = timetable_t::get_empty_value();//TODO:problem.get_empty_value()

	solution_t() {}

	solution_t(vector<int> &cities, unsigned price):path(cities), valid(true), cost(price) {
		//valid = true;
		//cost = price;
		//path = cities;
	}
};

//TODO: improve sparse data - find real flights first

//namespace greedy {
//solution_t greedy_search(problem_t &problem);
//}

void print_solution(const solution_t &solution, const problem_t &problem);
namespace dfs_improve {

bool in_vector(const vector<int> &vec, int elem) {
	//for (auto iter=vec.begin(); iter != vec.end(); ++iter) {
	for (unsigned i=0; i < vec.size(); ++i) {
		if (vec[i] == elem) {
			return true;
		}
	}
	return false;
}

solution_t recursive(const problem_t &problem, vector<int> &path, int day, unordered_set<int> cities_to_visit)
{
	solution_t best_solution;
	if (cities_to_visit.size() == 0) { //last flight, end recursion
		return solution_t(path, get_path_cost(path, problem));
	}
	for (auto i=cities_to_visit.cbegin(); i!=cities_to_visit.cend(); ++i) {
		vector<int> path_to_test = path;
		
		unordered_set<int> new_cities_to_visit = cities_to_visit;
		new_cities_to_visit.erase(*i);
		path_to_test[day] = *i;
		solution_t proposal=recursive(problem, path_to_test, day+1, new_cities_to_visit);
		if (proposal.cost < best_solution.cost) {
			best_solution = proposal;
		}
	}
	return best_solution;
}

solution_t dfs_improve(const solution_t &initial_solution, const problem_t &problem)
{
	const int MAX_SUBPATH_LENGTH = 8;
	//int i;
	const int problem_size = problem.problem_size();
	unordered_set<int> cities_to_visit; //cities_to_visit.reserve(problem_size+1);
	//int current_city = problem.starting_city_index;
	vector<int> path = initial_solution.path;//path.reserve(problem_size+1);
	
	for (int i=0; i<problem_size-2; i+=MAX_SUBPATH_LENGTH-1) {
		//vector<int> subpath;
		int subpath_length;
		if (problem_size - i < MAX_SUBPATH_LENGTH) {
			subpath_length = problem_size - i;
		} else {
			subpath_length = MAX_SUBPATH_LENGTH;
		}
		cities_to_visit.clear();
		for (int j=i+1; j<i+subpath_length; ++j) {
			cities_to_visit.insert(path[j]);
		}
		//un_set_remove(cities_to_visit, problem.starting_city_index);
		//subpath.push_back(initial_solution.path[i]);
		//vector<int> opt_subpath = recursive(problem, subpath, i+subpath_length, i);
		path = recursive(problem, path, i+1, cities_to_visit).path;
		//for (j=i+1; j<i+subpath_length; ++j) {
		//	path[j] = opt_subpath[j];
		//}
	}
	return solution_t(path, get_path_cost(path, problem));
}

}//namespace dfs


namespace improve {

unsigned get_subcost_2swap (const vector<int> &path, int from, int place1, int place2, int to, int day, const problem_t &problem) {
	const unsigned empty = problem.get_empty_value();
	if (problem.get_price(path[from], path[place1], day-1) == empty) {
		return empty;
	}
	if (problem.get_price(path[place1], path[place2], day) == empty) {
		return empty;
	}
	if (problem.get_price(path[place2], path[to], day+1) == empty) {
		return empty;
	}
	return problem.get_price(path[from], path[place1], day-1) + problem.get_price(path[place1], path[place2], day) + problem.get_price(path[place2], path[to], day+1);
}


solution_t swap_improve(const solution_t &initial_solution, const problem_t &problem) {
	unsigned i;
	const unsigned size = problem.problem_size();
	//vector<int> path = greedy::greedy_search(problem).path;
	/*vector<int> path(size+1);
	path[0] = path[size] = problem.starting_city_index;
	int city_to_add=0;
	for (i=1; i<size; ++i) {
		if (city_to_add == problem.starting_city_index) {
			++city_to_add;
		}
		path[i] = city_to_add;
		++city_to_add;
	}*/
	solution_t temp_sol = initial_solution;
	solution_t best_sol = temp_sol;
	if (size < 3) {
		return temp_sol;
	}

	bool found = true;
	while (found) {
		//cout << "testujem\n";
		found = false;
		for (i=1; i<size-1; ++i) {
			unsigned cur_subcost = get_subcost_2swap(temp_sol.path, i-1, i, i+1, i+2, i, problem);
			unsigned new_subcost = get_subcost_2swap(temp_sol.path, i-1, i+1, i, i+2, i, problem);
			//cout << "subtest cur:" << cur_subcost << " new:" << new_subcost << endl;
			if (new_subcost < cur_subcost) {
				//cout << "menim\n";
				found = true;
				int tmp = temp_sol.path[i+1];
				temp_sol.path[i+1] = temp_sol.path[i];
				temp_sol.path[i] = tmp;
				if (temp_sol.cost == problem.get_empty_value()) {
					//cout << "neplatna\n";
					temp_sol.cost = get_path_cost(temp_sol.path, problem);
				} else {
					//cout << "platna\n";
					temp_sol.cost -= (cur_subcost - new_subcost);
				}
				best_sol = temp_sol;
#ifdef DEBUG
				if (get_path_cost(temp_sol.path, problem) != temp_sol.cost) {
					cout << "BUG: scratch:" << get_path_cost(temp_sol.path, problem) << " diff: " << temp_sol.cost << "\n";
					exit(-1);
				}
#endif
			}
		}
	}

	return best_sol;
}

}//end swap
	
namespace anneal {

unsigned get_subcost_2swap (const vector<int> &path, int from, int place1, int place2, int to, int day, const problem_t &problem) {
	const unsigned empty = problem.get_empty_value();
	if (problem.get_price(path[from], path[place1], day-1) == empty) {
		return empty;
	}
	if (problem.get_price(path[place1], path[place2], day) == empty) {
		return empty;
	}
	if (problem.get_price(path[place2], path[to], day+1) == empty) {
		return empty;
	}
	return problem.get_price(path[from], path[place1], day-1) + problem.get_price(path[place1], path[place2], day) + problem.get_price(path[place2], path[to], day+1);
}


class annealing_schedule {
private:
	const struct timespec &start_time;
	struct timespec real_start;
	struct timespec planned_end;
	long elapsed_time;
	long real_interval;
	//const int interval;
	minstd_rand &generator;
	const problem_t &problem;
	//annealing_schedule(){}
	const double init_temp;
	double temperature;
	//uniform_real_distribution<double> dist(0.0, 1.0);
	uniform_real_distribution<double> dist;
	int skip_timer;

	void adjust_temp() {
		++skip_timer;
		if (skip_timer < 10000) {
			return;
		}
		skip_timer = 0;
		struct timespec now;
		clock_gettime(CLOCK_MONOTONIC_RAW, &now);
		elapsed_time = (now.tv_sec - real_start.tv_sec) * 1000*1000*1000 + now.tv_nsec - real_start.tv_nsec;
		temperature = (1 - (double(elapsed_time) / real_interval)) * init_temp;
		if (temperature < 0 /*|| elapsed_time > real_interval*/) {
			temperature = 0;
		}
		//cout << "elapsed: " << elapsed_time << " real int:" << real_interval << " elapsed/int: " << (double(elapsed_time) / real_interval)<< " temp: " << temperature << endl;
	}

public:
	annealing_schedule(const struct timespec &_start_time, const long _interval, minstd_rand &_generator, const problem_t &_problem, const double initial_temperature):
		start_time(_start_time), elapsed_time(0), /*interval(_interval),*/ generator(_generator), problem(_problem), init_temp(initial_temperature), temperature(initial_temperature), dist(uniform_real_distribution<double> (0.0,1.0)), skip_timer(0)  {
		if (clock_gettime(CLOCK_MONOTONIC_RAW, &real_start) != 0) {
			//TODO: what if clock fails?
#ifdef DEBUG
			cout << "CLOCK FAILED\n";
			exit (-1);
#endif
			real_start = start_time;
		}
		planned_end.tv_sec = start_time.tv_sec + _interval/(1000l*1000l*1000l);
		planned_end.tv_nsec = start_time.tv_nsec + _interval%(1000l*1000l*1000l);
		while (planned_end.tv_nsec > (1000l*1000l*1000l)) {
			++planned_end.tv_sec;
			planned_end.tv_nsec -= (1000l*1000l*1000l);
		}
		real_interval = (planned_end.tv_sec - real_start.tv_sec) * 1000l*1000l*1000l + planned_end.tv_nsec - real_start.tv_nsec; 
#ifdef DEBUG
		cout << "starting temp: " << init_temp << endl;
#endif
	}

	inline bool should_swap(const unsigned new_subcost, const unsigned cur_subcost) {
		adjust_temp();
		if (cur_subcost == problem.get_empty_value()) {
			//cout << "nedef\n";
			return true;
		}
		if (cur_subcost > new_subcost) {
			//cout << "stary " << cur_subcost << " vacsi nez " << new_subcost << endl;
			//cout << "lepsi\n";
			return true;
		}
		if (temperature == 0) {
			//cout << "konecna\n";
			return false;
		}
		if (exp(-( (double)new_subcost - cur_subcost) /temperature) > dist(generator)) {
			//cout << "horsi stastie, vzorcek: "<<exp(-( (double)new_subcost - cur_subcost) /temperature)<< "\n";
			return true;
		}
		//cout << "zly\n";
		return false;
	}

	inline bool should_continue() const {
		if (temperature <= 0) {
			return false;
		}
		return true;
	}
};

vector<int> ordered_initial_path(const problem_t &problem) {
	const unsigned size = problem.problem_size();
	vector<int> path(size+1);
	path[0] = path[size] = problem.starting_city_index;
	int city_to_add=0;
	for (unsigned i=1; i<size; ++i) {
		if (city_to_add == problem.starting_city_index) {
			++city_to_add;
		}
		path[i] = city_to_add;
		++city_to_add;
	}
	return path;
}

//TODO: nonlinear temperature, biggest discount search, piecewise brute force
solution_t anneal_search(const vector<int> &initial_path, const problem_t &problem, const int initial_temperature) {
	unsigned i;
	const unsigned size = problem.problem_size();
	minstd_rand generator;
	generator.seed(1);
	uniform_int_distribution<int> distribution(1, size-2);
	annealing_schedule  schedule(start_time, TIME_NS-TIME_NS_RESERVE, generator, problem, initial_temperature);
	//vector<int> path = greedy::greedy_search(problem).path;
	vector<int> path = initial_path;
	//vector<int> path(size+1);
	/*path[0] = path[size] = problem.starting_city_index;
	int city_to_add=0;
	for (i=1; i<size; ++i) {
		if (city_to_add == problem.starting_city_index) {
			++city_to_add;
		}
		path[i] = city_to_add;
		++city_to_add;
	}*/
	solution_t temp_sol(path, get_path_cost(path, problem));
	solution_t best_sol = temp_sol;
	if (size < 3) {
		return temp_sol;
	}
#ifdef DEBUG
	unsigned long dummy_counter = 0;
#endif
	//bool cont = true;
	while (schedule.should_continue()) {
#ifdef DEBUG
		++dummy_counter;
		if (dummy_counter == ULONG_MAX) {
			cout << "Counter overflow\n";
			dummy_counter = 0;
		}
#endif
		// cont = false;
		//cout << "testujem\n";
		//found = false;
		i = distribution(generator);
		//for (i=1; i<size-1; ++i) {
			unsigned cur_subcost = get_subcost_2swap(temp_sol.path, i-1, i, i+1, i+2, i, problem);
			unsigned new_subcost = get_subcost_2swap(temp_sol.path, i-1, i+1, i, i+2, i, problem);
			//cout << "subtest cur:" << cur_subcost << " new:" << new_subcost << endl;
			if (schedule.should_swap(new_subcost, cur_subcost)) {
				//found = true;
				int tmp = temp_sol.path[i+1];
				temp_sol.path[i+1] = temp_sol.path[i];
				temp_sol.path[i] = tmp;
				if (temp_sol.cost == problem.get_empty_value() || new_subcost == problem.get_empty_value()) {
					//cout << "neplatna\n";
					temp_sol.cost = get_path_cost(temp_sol.path, problem);
				} else {
					//cout << "platna\n";
					temp_sol.cost -= (cur_subcost - new_subcost);
				}
				if (best_sol.cost > temp_sol.cost) {
					best_sol = temp_sol;
				}
				//cout << temp_sol.cost << " menim\n";
#ifdef DEBUG
				if (get_path_cost(temp_sol.path, problem) != temp_sol.cost) {
					cout << "BUG: scratch:" << get_path_cost(temp_sol.path, problem) << " diff: " << temp_sol.cost << "\n";
					exit(-1);
				}
#endif
			}
			else {
			//	cout << "nemenim\n";
			}
		//}
	
	}
#ifdef DEBUG
	cout << "COunter: " << dummy_counter << " / " << ULONG_MAX << "  starting temperature: " << initial_temperature << " end state: " << temp_sol.cost << endl;
#endif
	return best_sol;
}

}//end anneal
	
namespace concorde {

solution_t process_concorde_output(const string &filename, const problem_t &problem) {
	ifstream in(filename);
	int size=0;
	in >> size;
	vector<int> raw_path(size);
	vector<int> path(size+1);
	if (size != problem.problem_size()) {
#ifdef DEBUG
		cout << "ERROR: concorde size mismatch\n";
		exit(-1);
#endif
		return solution_t();
	}
	for (int i=0; i<size; ++i) {
		in >> raw_path[i];
	}
	auto start_pos = find(raw_path.begin(), raw_path.end(), problem.starting_city_index);
#ifdef DEBUG
	if (start_pos == raw_path.end()) {
		cout << "ERROR: starting city not found\n";
		exit(-1);
	}
#endif
	for (int i=0; i<=size; ++i) {
		if (start_pos == raw_path.end()) {
			start_pos = raw_path.begin();
		}
		path[i] = *start_pos;
		++start_pos;
	}
	vector<int> reverse_path;
	reverse_path.reserve(path.size());
	for (auto iter=path.crbegin(); iter != path.crend(); ++iter) {
		reverse_path.push_back(*iter);
	}
	solution_t reverse_solution = improve::swap_improve(solution_t(reverse_path, get_path_cost(reverse_path, problem)), problem);
	solution_t orig_solution = improve::swap_improve(solution_t(path, get_path_cost(path, problem)), problem);
#ifdef DEBUG
	cout << "orig:" << orig_solution.cost << " reverse:" << reverse_solution.cost << "\n";
#endif

	if (orig_solution.cost > reverse_solution.cost) {
#ifdef DEBUG
		cout << "reverse is better\n";
#endif
		return reverse_solution;
	}
#ifdef DEBUG
		cout << "orig is better\n";
#endif
return orig_solution;
}

typedef enum {
	average,
	median
} averaging_method_t;

solution_t concorde_search(const averaging_method_t method, const problem_t &problem) {
	const unsigned size = problem.problem_size();
	const string filename = "k.tsp";//increase lenght with caution
	const string out_filename = "k.sol";//increase lenght with caution
	ofstream of(filename, ios_base::out|ios_base::trunc);
	of << "NAME: kiwi\nTYPE: TSP\nDIMENSION: ";
	of << size;
	of << "\nEDGE_WEIGHT_TYPE: EXPLICIT\nEDGE_WEIGHT_FORMAT: UPPER_ROW\nEDGE_WEIGHT_SECTION\n";

	if (method == average) {
		for (unsigned i=0; i<size-1; ++i) {
			for (unsigned j=i+1; j<size; ++j) {
				int count = 0;
				int sum = 0;
				for (unsigned k=0; k<size; ++k) {
					unsigned price = (problem.get_price(i, j, k));
					if (price != problem.get_empty_value()) {
						++count;
						sum += price;
					}
					price = (problem.get_price(j, i, k));
					if (price != problem.get_empty_value()) {
						++count;
						sum += price;
					}
				}
				of << sum / count << " ";
			}
			of << endl;
		}
	} else if (method == median) {
		for (unsigned city1=0; city1<size-1; ++city1) {
			for (unsigned city2=city1+1; city2<size; ++city2) {
				//weird method to compute median?
				priority_queue<int> sort;
				for (unsigned day=0; day<size; ++day) {
					if (problem.get_price(city1, city2, day) != problem.get_empty_value()) {
						sort.emplace(problem.get_price(city1, city2, day));
					}
					if (problem.get_price(city2, city1, day) != problem.get_empty_value()) {
						sort.emplace(problem.get_price(city2, city1, day));
					}
				}
				const unsigned sort_size = sort.size();
				for (unsigned i=0; i<sort_size/2;++i) {
					sort.pop();
				}
				of << sort.top() << " ";
			}
			of << endl;
		}
	}

	of << "EOF";
	of.close();

	//char *const argv[] = {"solv", "-x", "-s", "1", "kiwi.tsp", NULL};
	//char const *argvorig[] = {"solv", "-x", "-s", "1", "k.tsp", NULL};
	//TODO: something more elegant to prepare argv?
	const string executable = "./solv";
	const vector<string> argvorig =  {executable, "-x", "-s", "1", "-o", out_filename, filename};
	char *argvcopy[argvorig.size()+1];
	unsigned i=0;
	while (i < argvorig.size()) {
		argvcopy[i] = strdup(argvorig[i].c_str()); //TODO: memory leak?
		++i;	
	}
	argvcopy[argvorig.size()] = NULL;
	//char *const argv[] = {(argvcopy[0].c_str()),(argvcopy[1].c_str()),(argvcopy[2].c_str()),(argvcopy[3].c_str()),(argvcopy[4].c_str()),NULL};
	char *const envp[] = {NULL};
	int wstatus;
	pid_t pid = fork();
	if (pid == 0) {
#ifndef DEBUG
		close(STDOUT_FILENO);
		close(STDERR_FILENO);
		close(STDIN_FILENO);
#endif
		execve(executable.c_str(), argvcopy, envp);
#ifdef DEBUG
		cout << "EXECVE failed\n";
#endif
		exit(-1);
	} else if (pid > 1) {
		waitpid(pid, &wstatus, 0);
		if (WIFEXITED(wstatus)) {
			//cout << "OK\n";
			return process_concorde_output(out_filename, problem);
		} else {
#ifdef DEBUG
			cout << "CONCORDE not wifexited\n";
#endif
			return solution_t();
		}
		//dead return
		//return solution_t();
	}
	else {
#ifdef DEBUG
		cout << "FORK failed\n";
		return solution_t();
#endif
	}
	//execvp(executable.c_str(), argvcopy);
	//dead code here
	return solution_t();
}

}//end concorde

namespace dfs {

bool in_vector(const vector<int> &vec, int elem) {
	//for (auto iter=vec.begin(); iter != vec.end(); ++iter) {
	for (unsigned i=0; i < vec.size(); ++i) {
		if (vec[i] == elem) {
			return true;
		}
	}
	return false;
}

solution_t recursive(const problem_t &problem, vector<int> &path)
{
	solution_t best_solution;

	if ((signed)path.size() == problem.problem_size()) { //last flight, end recursion
		path.push_back(problem.starting_city_index);
		return solution_t(path, get_path_cost(path, problem));
	}
	for (int i=0; i<problem.problem_size(); ++i) {
		if (in_vector(path, i)) {
			continue;
		}
		vector<int> path_to_test = path;
		path_to_test.push_back(i);
		solution_t proposal=recursive(problem, path_to_test);
		if (proposal.cost < best_solution.cost) {
			best_solution = proposal;
		}
	}
	return best_solution;
}

solution_t dfs_search(const problem_t &problem)
{
	int i;
	const int problem_size = problem.problem_size();
	unordered_set<int> cities_to_visit; cities_to_visit.reserve(problem_size+1);
	int current_city = problem.starting_city_index;
	vector<int> path;//path.reserve(problem_size+1);
	
	for (i=0; i<problem_size; ++i) {
		cities_to_visit.insert(i);
	}
	un_set_remove(cities_to_visit, problem.starting_city_index);
	path.push_back(current_city);

	return recursive(problem, path);
}

}//namespace dfs




namespace greedy {

solution_t greedy_search(const problem_t &problem)
{
	int i;
	const int problem_size = problem.problem_size();
	unordered_set<int> cities_to_visit; cities_to_visit.reserve(problem_size);
	int current_city = problem.starting_city_index;
	int today=0;
	//uint16_t best_price;
	int best_city;
	vector<int> path;path.reserve(problem_size);
	
	for (i=0; i<problem_size; ++i) {
		cities_to_visit.insert(i);
	}
	un_set_remove(cities_to_visit, problem.starting_city_index);
	path.push_back(current_city);

	while (!cities_to_visit.empty()) {
		unsigned best_price = problem.get_empty_value();
		//cout << "today:" << today <<  " in city:" << current_city << endl; 
		for (auto city=cities_to_visit.begin(); city != cities_to_visit.end(); ++city) {
			unsigned cur_price = problem.get_price(current_city, *city, today);
			//cout <<" testing city:" << *city << " price: " << cur_price << endl;
			if (cur_price < best_price) {
				//cout << "chosen\n";
				best_price = cur_price;
				best_city = *city;			
			}
		}
		if (best_price == problem.get_empty_value()) {
#ifdef DEBUG
			cerr << "ERROR: flight not found\n";
#endif
			return solution_t();
		}
		path.push_back(best_city);
		current_city = best_city;
		un_set_remove(cities_to_visit, best_city);
		++today;
	}
#ifdef DEBUG
	if (problem.get_price(current_city, problem.starting_city_index, today) ==  problem.get_empty_value()) {
			cerr << "ERROR: last flight missing\n";
	}
#endif
	path.push_back(problem.starting_city_index);
	return solution_t(path, get_path_cost(path, problem));
}

}//namespace greedy

void print_solution(const solution_t &solution, const problem_t &problem) {
	if (!solution.valid) {
#ifdef DEBUG
		cout<<"ERROR: invalid solution to print\n";
#endif
		return;
	}
	cout << solution.cost;

	//for (auto iter=solution.path.begin(); iter !=solution.path.end(); iter++) {
	for (int day=0; day<problem.problem_size(); ++day) {
		cout << endl;
		cout << problem.get_city_name(solution.path[day]) << " ";
		cout << problem.get_city_name(solution.path[day+1]) << " ";
		cout << day << " ";
		cout << problem.get_price(solution.path[day], solution.path[day+1], day);
	}
	
}

struct global_optimum_t {
	solution_t solution;
	pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
	bool initialized = false;
};


//another nasty global
global_optimum_t global_optimum;
const problem_t *global_problem;

bool register_result(const solution_t &solution) {
	bool rc = false;
	if ((pthread_mutex_lock(&(global_optimum.mutex))) != 0) {
#ifdef DEBUG
		cout << "CAN NOT LOCK\n";
		//cout << status << endl;
		exit(-1);
#endif
	}
	//crit section
		if (global_optimum.initialized == false || (solution.cost < global_optimum.solution.cost)) {
			global_optimum.solution = solution;
			global_optimum.initialized = true;
			rc = true;
		}
	//crit section end
	if (pthread_mutex_unlock(&(global_optimum.mutex)) != 0) {
#ifdef DEBUG
		cout << "CAN NOT UNLOCK\n";
		exit(-1);
#endif
	}
	return rc;

}

void run_search(void * func(void *), void *data);

void *dfs_thread_wrapper(void *start) {
	solution_t *data = (solution_t *)(start);
	solution_t copy = *data;
	solution_t solution = improve::swap_improve(dfs_improve::dfs_improve(copy, *global_problem), *global_problem);
	if (register_result(solution)) {
#ifdef DEBUG
		cout << "dfs is da best: " << solution.cost << "\n";
	} else {
		cout << "dfs  is NOT da best: " << solution.cost << "\n";
#endif
	}
	return NULL;
}


void *concorde_thread_wrapper(void *dummy)
{
	solution_t solution = concorde::concorde_search(concorde::average, *global_problem);
	if (register_result(solution)) {
#ifdef DEBUG
		cout << "concorde average is da best: " << solution.cost << "\n";;
	} else {
		cout << "concorde average is NOT da best: " << solution.cost << "\n";;
#endif
	}

	solution_t greedy_solution = improve::swap_improve(greedy::greedy_search(*global_problem), *global_problem);
	if (register_result(greedy_solution)) {
#ifdef DEBUG
		cout << "greedy is da best: " << greedy_solution.cost << "\n";;
	} else {
		cout << "greedy is NOT da best: " << greedy_solution.cost << "\n";;
#endif
	}
	if (solution.cost > greedy_solution.cost) {
		solution = greedy_solution;
	}


	solution_t solution2 = concorde::concorde_search(concorde::median, *global_problem);
	if (register_result(solution2)) {
#ifdef DEBUG
		cout << "concorde median is da best: " << solution2.cost << "\n";;
	} else {
		cout << "concorde median is NOT da best: " << solution2.cost << "\n";;
#endif
	}
//TODO: anneal best solution
	//anneal the result of concorde
	if (solution.cost > solution2.cost) {
		solution = solution2;
	}
	//solution = dfs_improve::dfs_improve(solution, *global_problem);
	run_search(dfs_thread_wrapper, &solution);
	
	solution = improve::swap_improve(anneal::anneal_search(solution.path, *global_problem, CA_TEMP), *global_problem);
	if (register_result(solution)) {
#ifdef DEBUG
		cout << "concorde anneal is da best: " << solution.cost << "\n";;
	} else {
		cout << "concorde anneal is NOT da best: " << solution.cost << "\n";;
#endif
	}


	return dummy;
}	

void *anneal_thread_wrapper(void *temperature)
{
	const double temp = (*(double *)temperature);
	solution_t solution = improve::swap_improve(anneal::anneal_search(anneal::ordered_initial_path(*global_problem), *global_problem, temp), *global_problem);
	if (register_result(solution)) {
#ifdef DEBUG
		cout << "anneal " << temp << " is da best: " << solution.cost << "\n";
	} else {
		cout << "anneal " << temp << " is NOT da best: " << solution.cost << "\n";
#endif
	}
	return NULL;
}	

void run_search(void * func(void *), void *data) {
	//int status;
	pthread_t thread;
	pthread_create(&thread, NULL, func, data);
}

void print_global_optimum(const problem_t &problem) {
	if (pthread_mutex_lock(&(global_optimum.mutex))) {
#ifdef DEBUG
		cout << "CAN NOT PRINT\n";
#endif
		return;
	}
	print_solution(global_optimum.solution, problem);
	if (pthread_mutex_unlock(&(global_optimum.mutex))) {
#ifdef DEBUG
		cout << "CAN NOT UNPRINT\n";
		exit(-1);
#endif
	}

}

bool timeout(const struct timespec &start, const struct timespec &end) {
	long elapsed_ns = 1000000000l*(end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec;
	if (elapsed_ns >= TIME_NS) {
#ifdef DEBUG
		cout << "TIME ELAPSED: " << elapsed_ns << endl;
#endif
		return true;
	}
	return false;
}

int main()
{
	asm volatile("": : :"memory");
#ifdef DEBUG
	cout << "START" << endl << flush;
	if ((sizeof(unsigned short) < sizeof(uint16_t) ) ||  (sizeof(uint16_t) >= sizeof(unsigned))) {
		cerr << "ERROR: SMALL INTEGER" << endl;
		return -1;
	}
	if (clock_status != 0) {
		cout << "NO CLOCK\n";
		return -1;
	}
	//cout << "empty val: " << problem_t::get_empty_value()<<endl;
#endif
	const problem_t problem;
	global_problem = &problem;
	if (global_problem->problem_size() <= 11) {
		print_solution(dfs::dfs_search(*global_problem), *global_problem);
	} else {
			//print_solution(improve::swap_improve(anneal::anneal_search(problem), problem), problem);
		double temperature[] = {INITIAL_TEMP};
		//double temperature[] = {};
		for (unsigned i=0; i<(sizeof(temperature)/sizeof(temperature[0])); ++i) {
			run_search(anneal_thread_wrapper, &(temperature[i]));
		}
		run_search(concorde_thread_wrapper, NULL);
	//print_solution(dfs::dfs_search(problem), problem);
	//cout << endl;
	//print_solution(improve::swap_improve(greedy::greedy_search(problem), problem), problem);
	//cout << endl;
	//print_solution(improve::swap_improve(concorde::concorde_search(problem), problem), problem);
	//cout << endl;
	//problem.print_table();
	//cout << "starting city: " << problem.starting_city_string << endl;
		struct timespec end_time;
		clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
		sleep(TIME_NS/1000000000l-3+start_time.tv_sec-end_time.tv_sec);
#ifdef DEBUG
		cout << "koncim dlhe cakanie\n";
#endif
		const struct timespec nanotimeout = {0, 1000000};
		struct timespec dummy;
		do {
			nanosleep(&nanotimeout, &dummy);
			clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
		} while (!timeout(start_time, end_time));
	
		print_global_optimum(*global_problem);
	}
#ifdef DEBUG
	cout << endl;
	cout << "END" << endl << flush;
#endif
	return 0;
}
