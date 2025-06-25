#include <iostream>
#include <fstream>
#include "fock/ParticleDrawing.h"
#include "fock/ImportanceTruncation.h"
#include "fock/Integrator.h"

std::vector<std::string> process_args(int argc, char** argv)
{
    std::vector<std::string> args;
    for (int i = 1; i < argc; ++i)
        args.push_back(argv[i]);
    return args;
}

#define CHECK_END(cond, msg) if (cond) { std::cerr << msg << std::endl; return 1; }

int main(int argc, char** argv)
{
    double gN = 0.0;
    uint npart = 6;
    int momentum = 0;
    double ener_cutoff = 150.0;
    uint ndraws = 8;
    std::string filename = "wave_function.dat";
    std::vector<std::string> args = process_args(argc, argv);
    for (auto it = args.begin(); it != args.end(); ++it) {
        if (*it == "-gN") {
            CHECK_END(++it == args.end(), "Error: -gN requires an argument");
            gN = std::stod(*it);
        } else if (*it == "-npart") {
            CHECK_END(++it == args.end(), "Error: -npart requires an argument");
            npart = std::stoi(*it);
        } else if (*it == "-momentum") {
            CHECK_END(++it == args.end(), "Error: -momentum requires an argument");
            momentum = std::stoi(*it);
        } else if (*it == "-ener_cutoff") {
            CHECK_END(++it == args.end(), "Error: -ener_cutoff requires an argument");
            ener_cutoff = std::stod(*it);
        } else if (*it == "-filename") {
            CHECK_END(++it == args.end(), "Error: -filename requires an argument");
            filename = *it;
        } else if (*it == "-ndraws") {
            CHECK_END(++it == args.end(), "Error: -ndraws requires an argument");
            ndraws = std::stoi(*it);
        } else {
            std::cout << "Usage: " << argv[0] << " [-gN <gN>] [-npart <npart>] [-momentum <momentum>] [-ener_cutoff <ener_cutoff>] [-filename <filename>] [-ndraws <ndraws>]" << std::endl;
            return 1;
        }
    }


    Mbs::TruncationParameters params;
    params.ener_cutoff = ener_cutoff;

    Mbs::integrator_setup(10);
    Mbs::Workspace ws = Mbs::importance_truncation_scheme(gN/npart, npart, momentum, params);
    Mbs::FockLinear evec = ws.eigenvector_combination(0);

    std::ofstream file(filename, std::ios::binary);
    srand(time(NULL));
    for (uint i = 0; i < ndraws; ++i) {
        std::cout << "### Drawing particle " << i << " ###" << std::endl;
        std::vector<Mbs::RPhi> positions = Mbs::draw_positions(evec, 1e-4);
        std::vector<Complex> wave_function = Mbs::fill_wave_function(positions, evec, 201, 1e-4);
        file.write(reinterpret_cast<const char*>(wave_function.data()), wave_function.size() * sizeof(Complex));
        file.flush();
    }
    file.close();
    Mbs::integrator_cleanup();

    return 0;
}
