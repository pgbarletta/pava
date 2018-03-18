#include "chemfiles.hpp"
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    std::cerr << "Usage: ./einvert in.pdb out.pdb" << '\n'
              << "Exiting now." << '\n';
    return 0;
  }

  chemfiles::Trajectory in_trj(argv[1]);
  auto frm_cnt = in_trj.nsteps();

  std::cout << "Frame count: " << frm_cnt << '\n';
  chemfiles::Trajectory out_trj(argv[2], 'w');
  for (size_t i = frm_cnt - 1; i > 0; --i) {
    const auto frm = in_trj.read_step(i);
    out_trj.write(frm);
  }
}
