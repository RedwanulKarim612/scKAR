#include <bits/stdc++.h>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

typedef long long ll;

// Function to read CSV files in the given format and merge into a count matrix
void merge_count_matrix(const string& input_path, const string& output_file, int cell_count) {
    map<string, vector<double>> kmer_map;
    set<string> all_kmers;
    vector<string> files;
    
    // Get the list of CSV files
    for (const auto& entry : fs::directory_iterator(input_path)) {
        string file = entry.path();
        if (file.substr(file.find_last_of(".") + 1) == "csv") {
            files.push_back(file);
        }
    }
    
    sort(files.begin(), files.end()); // Ensure files are read in order
    size_t file_count = files.size();
    
    // Read data from each file and populate kmer_map
    for (const auto& filename : files) {
        ifstream file_reader(filename);
        if (!file_reader.is_open()) {
            cerr << "Error opening file: " << filename << endl;
            continue;
        }
        
        string line;
        while (getline(file_reader, line)) {
            stringstream ss(line);
            string kmer;
            getline(ss, kmer, ',');
            all_kmers.insert(kmer);

            vector<double> counts(cell_count, 0.0);
            string index_str, count_str;
            while (getline(ss, index_str, ',') && getline(ss, count_str, ',')) {
                int file_index = stoi(index_str);
                double count = stod(count_str);
                counts[file_index] = count;
            }
            kmer_map[kmer] = counts;
        }
        file_reader.close();
    }

    // Write the merged count matrix
    ofstream output_writer(output_file);
    if (!output_writer.is_open()) {
        cerr << "Error opening output file: " << output_file << endl;
        return;
    }
    
    // Write data
    for (const auto& kmer : all_kmers) {
        output_writer << kmer;
        vector<double>& counts = kmer_map[kmer];
        for (double count : counts) {
            output_writer << "," << count;
        }
        output_writer << "\n";
    }
    
    output_writer.close();
}

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage: " << argv[0] << " <adj_path> <output_file> <cell_count>" << endl;
        return 1;
    }
    
    string adj_path = argv[1];
    string output_file = argv[2];
    int cell_count = stoi(argv[3]);
 
    merge_count_matrix(adj_path, output_file, cell_count);
    return 0;
}
