#include<bits/stdc++.h>  
#include <sys/stat.h>
#include <filesystem>
using namespace std;

typedef long long ll;

int output_count = 0;

ll start_time = time(0);
double MAX_RAM_USAGE = 8.0*1024.0*1024.0; // 10GB

void write_to_file(map<string, vector<pair<int, double>>> &kmer_map, string path){
    ofstream file_writer(path + "adj_" + to_string(output_count++) + ".csv", ios::app);
    stringstream ss;
    ll bufferSize = 0; // Track buffer size
    ll maxBufferSize = 4194304000; // 4GB for example, adjust based on your needs
    // string buffer = "";

    for (const auto &kmer_pair : kmer_map){
        stringstream lineBuffer; 
        lineBuffer << kmer_pair.first;
        for (const auto &count_pair : kmer_pair.second){
            lineBuffer << "," << count_pair.first << "," << count_pair.second; 
        }
        lineBuffer << "\n";

        string line = lineBuffer.str();
        ll lineSize = line.size();
        ss << line;
        bufferSize += lineSize; // Increment buffer size
        // buffer += line;
        lineBuffer.str(std::string());
        lineBuffer.clear();

        // Check if buffer exceeds maxBufferSize and write to file if it does
        if (bufferSize >= maxBufferSize){
            printf("Flushing, buffer size: %lld\n", bufferSize);
            file_writer.write(ss.str().c_str(), ss.str().size());
            file_writer.flush();
            bufferSize = 0;
            ss.clear(); 
            ss.str(std::string());
        }
    }

    // Write any remaining content in the buffer
    if (bufferSize > 0){
        printf("Final flush, buffer size: %lld\n", bufferSize);
        file_writer.write(ss.str().c_str(), bufferSize);
    }
    file_writer.close();
    ss.clear();
}

double get_file_size(string file_name){
    return std::filesystem::file_size(file_name)/1024.0;
}

int main(int argc, char** argv){
    // ios_base::sync_with_stdio(false);
    // cin.tie(NULL);
    string path = argv[1];
    cout << path << endl;
    ifstream file_mapping(path + string("/tpm_sum.csv"));
    map<string, int> file_map;
    string file_name;
    int file_no;
    // ll flush_count = stoll(argv[1]);
    ll file_no_count = 0;

    string adj_dir = path + string("/adj");
    mkdir(adj_dir.c_str(), 0777);

    while(true) {
        file_mapping >> file_name;
        if(file_mapping.eof()) {
            break;
        }
        file_mapping >> file_no;
        // cout << file_name << " " << file_no << "\n";
        file_map[file_name] = file_no_count;
        file_no_count++;
    }

    cout << file_map.size() << endl;

    map<string, vector<pair<int, double>>> kmer_map;
    int processed_count = 0;
    double cumulative_ram_usage = 0.0;
    for(auto it=file_map.begin(); it!=file_map.end(); it++){
        file_name = it->first + string("_1_filtered.csv");
        cout << file_name << endl;
        cout << path + string("/jellyfish/") + file_name << endl;
        ifstream file_reader(path + string("/jellyfish/") + file_name);
        string line;

        // if(processed_count < 96) {
        //     // cout << "Skipping: " << processed_count << endl;
        //     processed_count++;
        //     continue;
        // }

        while(getline(file_reader, line)){
            stringstream ss(line);
            string kmer, count;
            getline(ss, kmer, ',');
            getline(ss, count);
            kmer_map[kmer].push_back(make_pair(it->second, atof(count.c_str())));
        }
        cumulative_ram_usage += get_file_size(path + string("/jellyfish/") + file_name);
        cout << "Cumulative ram usage: " << cumulative_ram_usage << endl;
        // if(!remove((path + string("/jellyfish/") + file_name).c_str())) {
        //     cout << "File deleted successfully" << endl;
        // } else {
        //     cout << "Error deleting file" << endl;
        // }
        
        processed_count++;
        printf("Current: %d\n", processed_count);
        // if(processed_count % flush_count == 0){
        //     printf("file no: %d\n", processed_count);
        //     write_to_file(kmer_map, adj_dir + string("/"));
        //     kmer_map.clear();
        // }
        if (cumulative_ram_usage > MAX_RAM_USAGE){
            printf("file no: %d\n", processed_count);
            write_to_file(kmer_map, adj_dir + string("/"));
            kmer_map.clear();
            cumulative_ram_usage = 0.0;
        }
    }

    if (kmer_map.size() > 0){
        printf("Last file no: %d\n", processed_count);
        write_to_file(kmer_map, adj_dir + string("/"));
    }
    return 0;
}

// args: path_to_files