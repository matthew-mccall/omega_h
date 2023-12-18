#include <Omega_h_element.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_tag.hpp>
#include <Omega_h_array_ops.hpp>

template <typename T>
void printTagInfo(Omega_h::Mesh mesh, std::ostringstream& oss, int dim, int tag, std::string type) {
    auto tagbase = mesh.get_tag(dim, tag);
    auto array = Omega_h::as<T>(tagbase)->array();

    Omega_h::Real min = get_min(array);
    Omega_h::Real max = get_max(array);

    oss << "(" << tagbase->name().c_str() << ", " << dim << ", " << type << ")\n";
    oss << "\tNum Components: " << tagbase->ncomps() << "\n";
    oss << "\tMin, Max: " << min << ", " << max << "\n\n";
}

template <typename T>
void printValue(Omega_h::Mesh mesh, std::string tagname, int dim, int value) {
    auto array = mesh.get_array<T>(dim, tagname);
    auto each_eq_to = Omega_h::each_eq_to<T>(array, value);
    auto num = Omega_h::get_sum(each_eq_to);
    std::cout << "Num entities: " << num << "\n";
}

int main(int argc, char** argv)
{
    auto lib = Omega_h::Library(&argc, &argv);
    auto comm = lib.world();
    Omega_h::Mesh mesh = Omega_h::read_mesh_file(argv[1], lib.world());

    auto verbose = false;
    if (argc == 3) verbose = (std::string(argv[2]) == "on");

    const int rank = comm->rank();

    std::array<Omega_h::GO, 4> counts;
    for(int dim=0; dim < mesh.dim(); dim++)
       counts[dim] = mesh.nglobal_ents(dim);

    std::array<double, 4> imb;
    for(int dim=0; dim < mesh.dim(); dim++)
       imb[dim] = mesh.imbalance(dim);

    std::ostringstream oss;
    // always print two places to the right of the decimal
    // for floating point types (i.e., imbalance)
    oss.precision(2);
    oss << std::fixed;

    if(!rank) {
        oss << "\nMesh Entity Type: " << Omega_h::topological_singular_name(mesh.family(), mesh.dim()-1) << "\n";

        oss << "\nGlobal Mesh Entity Count and Imbalance (max/avg): (Dim, Entity Count, Imbalance)\n";
        for(int dim=0; dim < mesh.dim(); dim++)
            oss << "(" << dim << ", " << counts[dim] << ", " << imb[dim] << ")\n";

        oss << "\nTag Properties: (Name, Dim, Type)\n";
        for (int dim=0; dim < mesh.dim(); dim++)
        for (int tag=0; tag < mesh.ntags(dim); tag++) {
            auto tagbase = mesh.get_tag(dim, tag);
            if (tagbase->type() == OMEGA_H_I8)
                printTagInfo<Omega_h::I8>(mesh, oss, dim, tag, "I8");
            if (tagbase->type() == OMEGA_H_I32)
                printTagInfo<Omega_h::I32>(mesh, oss, dim, tag, "I32");
            if (tagbase->type() == OMEGA_H_I64)
                printTagInfo<Omega_h::I64>(mesh, oss, dim, tag, "I64");
            if (tagbase->type() == OMEGA_H_F64)
                printTagInfo<Omega_h::Real>(mesh, oss, dim, tag, "F64");
        }

        std::cout << oss.str();
    }

    if(verbose) {
        comm->barrier(); // write the per-part data at the end
        if(!rank) {
          std::cout << "\nPer Rank Mesh Entity Count: (Rank: Entity Count by Dim <0,1,2,3>)\n";
        }
        comm->barrier(); // write the per-part data at the end
        oss.str(""); // clear the stream

        std::array<Omega_h::LO, 4> counts = {0,0,0,0};
        for(int dim=0; dim < mesh.dim(); dim++)
            counts[dim] = mesh.nents(dim);
        oss << "(" << rank << ": " << counts[0] << ", "
                                        << counts[1] << ", "
                                        << counts[2] << ", "
                                        << counts[3] << ")\n";
        std::cout << oss.str();
    }

    std::string input;
    while(true) {
        std::cout << "Get num entities with value: options [tagname dim value] or [exit]\n";

        std::string name = "";
        int dim=0;
        int value=0;
        std::cin >> name;
        if (name == "exit") break;
        std::cin >> dim >> value;
        auto tagbase = mesh.get_tagbase(dim, name);
        if (tagbase->type() == OMEGA_H_I8)
            printValue<Omega_h::I8>(mesh, name, dim, value);
        if (tagbase->type() == OMEGA_H_I32)
            printValue<Omega_h::I32>(mesh, name, dim, value);
        if (tagbase->type() == OMEGA_H_I64)
            printValue<Omega_h::I64>(mesh, name, dim, value);
        if (tagbase->type() == OMEGA_H_F64)
            printValue<Omega_h::Real>(mesh, name, dim, value);
    }
    return 0;
}
