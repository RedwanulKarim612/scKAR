#include<bits/stdc++.h>
#include "./f_oneway_test.cpp"
#include "./bh_correction.cpp"

using namespace std;
namespace fs = std::filesystem;

// #define INF 99999999
typedef long long ll;

ll output_count = 1;


vector<string> string_separator(string str, string delimiter) {
    vector<string> tokens;
    stringstream check1(str);
    string intermediate;
    while(getline(check1, intermediate, ',')) {
        tokens.push_back(intermediate);
    }
    return tokens;
}

void pad_vector_with_zeros(vector<double> &vec, int count) {
    ll sz= vec.size();
    while((sz++) < count) {
        vec.push_back(0);
    }
}

void f_oneway_worker(
    const std::map<string, vector<pair<int, double>>>::iterator start, 
    const std::map<string, vector<pair<int, double>>>::iterator end, 
    vector<int> &mapping, 
    vector<ll> &cluster_cell_count, 
    ll total_cell_count,
    map<string, double> &p_values_mapping, // Pass the vector to collect p-values
    std::mutex &p_value_mutex, // Mutex for thread-safe access
    vector<ll> &tpm,
    ll prefilter_rowsum_threshold
) {
    cout << "Inside f_oneway_worker\n";
    // This function works on a subset of kmer_map from start to end
    for(auto it = start; it != end; ++it) {
        vector<vector<double>> groups = vector<vector<double>>(cluster_cell_count.size(), vector<double>());
        double prefilter_rowsum = 0;
        for(auto vec: it->second) {
            groups[mapping[vec.first]].push_back(vec.second);
            prefilter_rowsum += (ll)round(((tpm[vec.first] * vec.second) / 1000000))
        }
        // prefilter
        if(prefilter_rowsum < prefilter_rowsum_threshold) {
            it->second = vector<pair<int, double>>();
            p_value_mutex.lock();
            p_values_mapping[it->first] = INF;
            p_value_mutex.unlock();
        }
        else {
            // f_test
            for(int i=0;i<groups.size();i++) {
                pad_vector_with_zeros(groups[i], cluster_cell_count[i]);
            }
            AnovaResult res = f_oneway(groups);
            p_value_mutex.lock();
            p_values_mapping[it->first] = res.p_value;
            p_value_mutex.unlock();
        }
    }
    cout << "Exiting f_oneway_worker\n";
}

void f_oneway_on_map(map<string, vector<pair<int, double>>> &kmer_map,
    vector<int> &mapping,
    vector<ll> &cluster_cell_count,
    ll total_cell_count,
    vector<ll> &tpm,
    ll prefilter_rowsum_threshold
) {
    const int num_threads = 10;
    auto it = kmer_map.begin();
    std::vector<std::thread> threads;
    map<string, double> p_values_map; // Vector to collect p-values
    std::mutex p_value_mutex; // Mutex for thread-safe p-value collection

    // Calculate the number of kmers each thread should process
    size_t per_thread = kmer_map.size() / num_threads;
    
    for(int i = 0; i < num_threads; ++i) {
        auto start = it;
        std::advance(it, std::min(per_thread, static_cast<size_t>(std::distance(it, kmer_map.end()))));
        auto end = it;

        // Launch a new thread to process [start, end) interval of kmer_map
        threads.push_back(std::thread(f_oneway_worker, start, end, std::ref(mapping), std::ref(cluster_cell_count), total_cell_count, std::ref(p_values_map), std::ref(p_value_mutex), std::ref(tpm), prefilter_rowsum_threshold));
    }

    // Wait for all threads to complete
    for(auto &t : threads) {
        t.join();
    }
    printf("Starting Benjamini-Hochberg correction\n");
    vector<double> p_values_vector;
    for(auto ele: kmer_map) {
        if(p_values_map[ele.first] != INF) {
            // not prefiltered
            p_values_vector.push_back(p_values_map[ele.first]);
        }
    }
    vector<double> padj = compute_bh_correction(p_values_vector);
    ll i = 0;
    for(auto it = kmer_map.begin(); it != kmer_map.end(); ++it) {
        if(p_values_map[it->first] != INF) {
            // not prefiltered
            if(padj[i++] >= 0.05) {
                it->second = vector<pair<int, double>>();
            }
        }
        else {
            it->second = vector<pair<int, double>>();
        }
    }
    printf("Benjamini-Hochberg correction done\n");
}


void write_to_file_as_matrix(map<string, vector<pair<int, double>>> &kmer_map, string path, ll total_cell_count) {
    ofstream file_writer(path + "f_test_result_" + to_string(output_count++) + ".csv", ios::app | ios::out);
    stringstream ss;
    ll bufferSize = 0; // Track buffer size
    ll maxBufferSize = 4194304000; // 4GB for example, adjust based on your needs

    for (auto ele: kmer_map) {
        if(ele.second.size() == 0) {
            // didn't pass f test
            continue;
        }

        sort(ele.second.begin(), ele.second.end());
        stringstream lineBuffer;
        lineBuffer << ele.first;
        ll curr_file_no = 0;
        for(auto vec: ele.second) {
            while(curr_file_no < vec.first) {
                lineBuffer << ",";
                curr_file_no++;
            }
            lineBuffer << "," << vec.second;
            curr_file_no++;
        }
        while(curr_file_no < total_cell_count) {
            lineBuffer << ",";
            curr_file_no++;
        }
        lineBuffer << "\n";

        string line = lineBuffer.str();
        ll lineSize = line.size();
        ss << line;
        bufferSize += lineSize; 
        
        if (bufferSize >= maxBufferSize){
            printf("Flushing, buffer size: %lld\n", bufferSize);
            file_writer.write(ss.str().c_str(), bufferSize);
            bufferSize = 0; // Reset bufferSize
            ss.clear(); 
            ss.str(""); // Clear the contents of stringstream
        }
    }

    if (bufferSize > 0){
        printf("Final flush, buffer size: %lld\n", bufferSize);
        file_writer.write(ss.str().c_str(), bufferSize);
    }
    file_writer.close();
    ss.clear();
    ss.str("");
}

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


ll read_file_till_lineCount(map<string, vector<pair<int, double>>> &kmer_map, ll lineCount, ifstream &file_reader) {
    // cout << "Yooo inside!!\n";
    string line;
    ll curr_line = 0;
    ll bufferSize = 0;
    // file_reader.seekg(seekg_position);
    while (curr_line < lineCount && getline(file_reader, line)){
        bufferSize += line.size()+1;
        // cout << line << endl;
        stringstream ss(line);
        string kmer;
        getline(ss, kmer, ',');
        string file_no, count;
        while (getline(ss, file_no, ',') && getline(ss, count, ',')) {
            kmer_map[kmer].push_back({stoi(file_no), stod(count)});
        }
        // cout << count << endl;
        curr_line++;
    }
    // file_reader.close();
    return bufferSize;
}

ll read_file_till_kmer(map<string, vector<pair<int, double>>> &kmer_map, string till_kmer, string filename, ll seekg_position) {
    ifstream file_reader(filename);
    // cout << filename << endl; 
    string line;
    ll bufferSize = 0;
    file_reader.seekg(seekg_position);
    while (getline(file_reader, line)){
        bufferSize += line.size()+1;
        stringstream ss(line);
        string kmer;
        getline(ss, kmer, ',');
        if(kmer > till_kmer) {
            bufferSize -= line.size()+1;
            break;
        }
        string file_no, count;
        while (getline(ss, file_no, ',') && getline(ss, count, ',')) {
            kmer_map[kmer].push_back({stoi(file_no), stod(count)});
        }
    }
    file_reader.close();
    return bufferSize;
}

ll readLeidenMapping(string file, vector<ll> &tpm, vector<int> &mapping) {
    ifstream file_reader(file);
    string line;
    string cell_barcode, cell_tpm, cluster_id;
    while(getline(file_reader, line)) {
        stringstream ss(line);
        ss >> cell_barcode >> cell_tpm >> cluster_id;
        tpm.push_back(stoll(cell_tpm));
        mapping.push_back(stoi(cluster_id));
    }
    return tpm.size();
}


int main(int argc, char** argv) {
    string base_path = argv[1];
    string path = base_path + "/adj/";
    string output_path = base_path + "/f_test_results/";
    string leiden_file = base_path + "/cluster_tpm.csv";

    // if no dir output_path exists, create it
    if (!fs::exists(output_path)) {
        fs::create_directory(output_path);
    }

    ll lines_to_read = stoll(argv[2]);
    ll prefilter_rowsum_threshold = stoll(argv[3]);

    vector<ll> tpm = vector<ll>();
    vector<int> mapping = vector<int>();

    ll total_cell_count = readLeidenMapping(leiden_file, tpm, mapping);
    cout << "Total cell count: " << total_cell_count << endl;
    ll total_cluster_count = *max_element(mapping.begin(), mapping.end()) + 1;
    vector<ll> cluster_cell_count = vector<ll>(total_cluster_count, 0);
    for(auto ele: mapping) {
        cluster_cell_count[ele]++;
    }
    cout << "Total cluster count: " << total_cluster_count << endl;
    for(auto ele: cluster_cell_count) {
        cout << ele << " ";
    }
    cout << "---Starting Batch---\n";
    vector<string> files = get_list_files(path, "csv");

    string base_file = files[0];
    vector<ll> last_seekg_postions = vector<ll>(total_cell_count, 0);
    
    
    ifstream base_file_reader(base_file);

    map<string, vector<pair<int, double>>> total_map;

    // as long as there are content to read in base file
    int iter = 0;
    while(base_file_reader.peek() != EOF) {
        // ll base_seekg_pos = base_file_reader.tellg();
        read_file_till_lineCount(total_map, lines_to_read, base_file_reader);
        string till_kmer = total_map.rbegin()->first;
        for(int i=1;i<files.size();i++) {
            ll bytes_read = read_file_till_kmer(total_map, till_kmer, files[i], last_seekg_postions[i]);
            last_seekg_postions[i] += bytes_read;
        }
        // run f test
        cout << "Running f test\n";
        f_oneway_on_map(total_map, mapping, cluster_cell_count, total_cell_count, tpm, prefilter_rowsum_threshold);
        cout << "F test done\n";
        cout << "Writing to file\n";
        write_to_file_as_matrix(total_map, output_path, total_cell_count);
        cout << "Writing to file done\n";
        total_map.clear();
        cout << "Iter: " << iter++ << " " << output_count-1 <<endl;
    }

    string last_till_kmer = "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ";
    for(int i=1;i<files.size();i++) {
        read_file_till_kmer(total_map, last_till_kmer, files[i], last_seekg_postions[i]);
    }
    cout << "Running f test on last chunk\n";
    f_oneway_on_map(total_map, mapping, cluster_cell_count, total_cell_count, tpm, prefilter_rowsum_threshold);
    cout << "F test done\n";
    cout << "Writing to file in the last chunk\n";
    write_to_file_as_matrix(total_map, output_path, total_cell_count);

    base_file_reader.close();
    cout << "Lastly you see me!! YOOOO\n";
}



// args: args path, lines_to_read, prefilter_rowsum_threshold
// generally lines_to_read = 700000 for 4k cells in our local machine
// renal cell pre_filter = 20
// axolotl cell pre_filter = 100
