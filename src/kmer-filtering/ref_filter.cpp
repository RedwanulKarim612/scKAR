#include<bits/stdc++.h>
using namespace std;
namespace fs = std::filesystem;
typedef long long ll;
ll maxBufferSize = 4194304000/4; //1GB


vector<string> get_list_files(string path, string ext) {
    vector<string> files;
    for (const auto & entry : fs::directory_iterator(path)){
        string file = entry.path();
        // csv extension and starts with "adj"
        if ((file.substr(file.find_last_of(".") + 1) == ext) && (file.find("adj") != string::npos)){
            files.push_back(file);
        }
    }
    return files;
}

vector<string> string_separator(string str, string delimiter) {
    vector<string> tokens;
    stringstream check1(str);
    string intermediate;
    while(getline(check1, intermediate, ',')) {
        tokens.push_back(intermediate);
    }
    return tokens;
}

void write_to_file(stringstream& ss, string adj_file){
    string adj_file_name = adj_file.substr(0, adj_file.find_last_of("."));
    adj_file_name += ".filtered.csv";
    ofstream out_file(adj_file_name, ios::app | ios::out);
    ll buffer_size = ss.str().size();
    out_file.write(ss.str().c_str(), buffer_size);
    // out_file.flush();
    ss.str("");
    ss.clear();
    out_file.close();
}

void ref_filtering_worker(vector<string>& ref_data, string adj_file){
    cout << "entering ref_filtering_worker for " << adj_file << "\n";
    ifstream adj_file_stream(adj_file);
    string ref_line;
    string adj_line;
    vector<string> adj_data;
    int ref_index = 0;
    stringstream ss;
    ll ref_data_size = ref_data.size();
    ll cur_size = 0;

    while(getline(adj_file_stream, adj_line)){
        vector<string> adj_tokens = string_separator(adj_line, ",");
        string adj_kmer = adj_tokens[0];
        string ref_kmer = ref_data[ref_index];
        if (ref_index >= ref_data_size){
            ss << adj_line << '\n';     // all ref_kmers are exhausted
            cur_size += adj_line.size();
        }
        else if (adj_kmer == ref_kmer){
            ref_index++;   //adj_kmer in reference so filtering out
        }
        else if (adj_kmer < ref_kmer){
            // add adj_line to ss
            ss << adj_line << '\n';   //adj_kmer not in reference so keeping
            cur_size += adj_line.size();
        }
        else{
            while(ref_index < ref_data_size && adj_kmer > ref_data[ref_index]){
                ref_index++;
            }
            ref_kmer = ref_data[ref_index];
            if (adj_kmer == ref_kmer){
                ref_index++;
            }
            else{
                ss << adj_line << '\n';   //adj_kmer not in reference so keeping
                cur_size += adj_line.size();
            }
        }
        if(cur_size >= maxBufferSize){
            cout << "writing to file: " << adj_file << "\n"; 
            write_to_file(ss, adj_file);
            cur_size = 0;
        }
    }
    if(cur_size > 0){
        cout << "Final writing to file: " << adj_file << "\n";
        write_to_file(ss, adj_file);
    }
    cout << "exiting ref_filtering_worker for: " << adj_file << "\n";

}

void ref_filtering_manager(string ref_file, vector<string> adj_files, long long num_threads ){
    vector<string> ref_data;
    ifstream ref_file_stream(ref_file);
    string ref_line;
    while(getline(ref_file_stream, ref_line)){
        ref_data.push_back(ref_line);
    }

    std::vector<std::thread> threads;
    int thread_count = 0;
    for(int i=0;i<adj_files.size();i++){
        threads.push_back(std::thread(ref_filtering_worker, std::ref(ref_data), adj_files[i]));
        thread_count++;

        if(thread_count == num_threads){
            for(auto &t : threads){
                t.join();
            }
            threads.clear();
            thread_count = 0;
        }
    }
    for(auto &t : threads){
        t.join();
    }
    threads.clear();
    cout << "exiting ref_filtering_manager\n";
}


int main(int argc, char* argv[]){
    string folder_path = argv[1];
    string ref_file = argv[2];
    long long num_threads = stoll(argv[3]);

    ifstream ref_file_stream(ref_file);
    vector<string> adj_files = get_list_files(folder_path + "/adj/", "csv");


    cout << adj_files.size() << endl;
    for(int i=0;i<adj_files.size();i++){
        cout << adj_files[i] << endl;
    }

    auto start = chrono::high_resolution_clock::now();
    ref_filtering_manager(ref_file, adj_files, num_threads);
    auto end = chrono::high_resolution_clock::now();
    cout << "time taken: " << chrono::duration_cast<chrono::minutes>(end - start).count() << " minutes" << endl;

    return 0;
}

// args: folder_path, ref_T_file, num_threads