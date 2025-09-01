# AGENT INSTRUCTIONS

## Code style
- Format all C++ sources with `clang-format -style=llvm -i`.
- Always add comprehensive comments for documentation and maintainability:
  - **Function-level comments**: Document the purpose, parameters, return values, and any important assumptions or constraints for every function
  - **Complex logic comments**: Add detailed explanations for non-trivial algorithms, mathematical computations, regex patterns, and scientific formulations
  - **Code section comments**: Use section headers (e.g., `// ========== VALIDATION ==========`) to organize large functions into logical blocks
  - **Inline comments**: Explain why specific decisions were made, especially for scientific or numerical considerations
  - **Example transformations**: Provide concrete examples in comments when dealing with data transformations or chemical reactions

## Dependencies
- Ensure the Eigen3 library is installed before building:
  `apt-get update && apt-get install -y libeigen3-dev`

## Testing
Before committing changes run:
- `cmake -S . -B build`
- `cmake --build build`
- `ctest --test-dir build`
