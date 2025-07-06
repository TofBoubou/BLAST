// Add this debug output in coefficient_calculator.cpp line 113-117:
        std::vector<double> c_local(n_species);
        for (std::size_t j = 0; j < n_species; ++j) {
            c_local[j] = inputs.c(j, i);
            std::cout << "eta=" << i << " species=" << j << " c=" << c_local[j] << std::endl;
        }
EOF < /dev/null
