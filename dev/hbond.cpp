#include <chemfiles.hpp>
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: ./ehbond in.pdb out" << '\n' << "Exiting now." << '\n';
    return 0;
  }
  chemfiles::Trajectory in_trj(argv[1]);
  auto frm_cnt = in_trj.nsteps();
  std::cout << "Frame count: " << frm_cnt << '\n';

  for (size_t i = 0; i < frm_cnt; ++i) {
    const auto frm = in_trj.read_step(i);
    const chemfiles::Topology top = frm.topology();
    const auto xyz = frm.positions();

    const auto natoms = top.size();
    for (size_t j = 0; j < natoms - 1; j++) {
      for (size_t k = j; k < natoms; k++) {
        std::cout << "  " << frm.distance(j, k) << '\n';
      }
    }
  }

  return 0;
}
