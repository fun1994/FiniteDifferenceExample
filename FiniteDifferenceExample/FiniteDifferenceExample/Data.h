#pragma once
#include <vector>
#include <fstream>

class Data {
	void save(std::vector<double>& data, std::string path, std::string filename) {
		std::ofstream file("./data/" + path + "/" + filename + ".txt");
		for (int i = 0; i < data.size(); i++) {
			file << data[i];
			if (i < data.size() - 1) {
				file << " ";
			}
		}
		file.close();
	}
public:
	std::vector<double> x_exact;
	std::vector<double> phi_exact;
	std::vector<double> x_calculated;
	std::vector<double> phi_calculated;
	void save(std::string path, std::string convection) {
		save(x_exact, path, "x_exact");
		save(phi_exact, path, "phi_exact");
		save(x_calculated, path, "x_calculated_" + convection);
		save(phi_calculated, path, "phi_calculated_" + convection);
	}
};
