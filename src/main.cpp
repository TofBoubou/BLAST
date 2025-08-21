#include "blast/core/application_runner.hpp"

int main(int argc, char* argv[]) {
  blast::core::ApplicationRunner app;
  auto result = app.run(argc, argv);
  return result.exit_code;
}