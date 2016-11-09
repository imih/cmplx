#ifndef MPI_COMMON_H
#define MPI_COMMON_H

namespace cmplx {
enum MessageType {
  SIMUL_PREREQUEST,
  SIMUL_REQUEST,
  SIMUL_RESPONSE,
  SIMUL_END,
  SIMUL_PARAMS
};
}

#endif  // MPI_COMMON_H
