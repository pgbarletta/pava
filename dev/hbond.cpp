#include <chemfiles.hpp>
#include <iostream>
#include <string>
#include <unordered_map>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: ./ehbond in.pdb out" << '\n' << "Exiting now." << '\n';
    return 0;
  }
  chemfiles::Trajectory in_trj(argv[1]);
  auto frm_cnt = in_trj.nsteps();
  std::cout << "Frame count: " << frm_cnt << '\n';

  // for (size_t i = 0; i < frm_cnt; ++i) {
  //   const auto frm = in_trj.read_step(i);
  //   const chemfiles::Topology top = frm.topology();
  //   const auto xyz = frm.positions();
  //
  //   const auto natoms = top.size();
  //   for (size_t j = 0; j < natoms - 1; j++) {
  //     for (size_t k = j; k < natoms; k++) {
  //       std::cout << "  " << frm.distance(j, k) << '\n';
  //     }
  //   }
  // }

  std::unordered_map<std::string, int> mapa = {
      {"uno", 1}, {"dos", 2}, {"tre", 3}};

  for (std::unordered_map<std::string, int>::const_iterator ite = mapa.cbegin();
       ite != mapa.cend(); ++ite) {
    std::cout << ite->first << "  " << ite->second << '\n';
  }

  chemfiles::Atom atomo("NOS", "N");
  atomo.set("masa", 14);
  atomo.set("largo", 13);
  atomo.set("ancho", 12);

  std::cout << atomo.get("masa")->as_double() << "  " << '\n';
  std::cout << atomo.get("largo")->as_double() << "  " << '\n';
  std::cout << atomo.get("ancho")->as_double() << "  " << '\n';

  auto ite = atomo.properties_begin();
  auto ult = atomo.properties_end();

  for (ite; ite != ult; ++ite) {
    std::cout << ite->first << "  " << ite->second.as_double() << '\n';
  }

  const int count = std::stoi(argv[2]);
  std::cout << "count:  " << count << '\n';

  for (size_t i = 0; i < count; i++) {
    std::cout << atomo.get_each().first << "  "
              << atomo.get_each().second.as_double() << '\n';
  }

  return 0;
}
