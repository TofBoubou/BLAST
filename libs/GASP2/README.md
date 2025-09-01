# GASP² : Gas–Ablation & Surface / Pyrolysis Phenomena

GASP² is an open-source framework for modeling gas–surface and bulk interactions in chemically reacting high-enthalpy flows. The library presently targets catalysis calculations driven by an XML description of the surface chemistry.

## Building and Testing

The project uses CMake and requires a C++23 capable compiler.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

Running the tests is recommended to verify that the build and example input files behave as expected.

## Gamma Model

### Input File

To evaluate surface chemistry with recombination probabilities, set `<reactions model="gamma_model">` inside the XML file placed in the `inputs` directory. Important elements:

1. `<gasp2 surface_model="catalysis">` root element (use `non_catalytic` to disable reactions).
2. Optional `<first_order>true</first_order>` enabling a simplified first-order theory (required for `GammaConsistent` reactions).
3. A sequence of `<reaction>` entries. Each `reaction` has a `type` attribute selecting `GammaGiven`, `GammaT`, or `GammaConsistent`.
4. Within a reaction, supply either `<gammas>` or `<gammaT>` coefficients:
   * `<gammas>` entries like `Species:Value`, separated by commas or semicolons; values must lie in (0,1).
   * `<gammaT>` entries like `Species:A:E:Tmin:Tmax` giving an Arrhenius expression; `Tmax = -1` removes the upper validity bound.

### Supported Reactions and Parameters

The `gamma_model` covers heterogeneous recombination of at most two reactants yielding one product. Reactions are mass balanced automatically. Validation enforces:

* Every reactant listed in the `formula` has an associated gamma/gammaT entry and no coefficients are supplied for non-reactants.
* Each reaction uses either constant or temperature-dependent coefficients exclusively.
* After normalization duplicates are rejected and per-species gammas across reactions may not exceed unity.
* For the special `Rini_model`, a global `<gamma_w>` parameter in `[0,1]` is required.

### Initialize

```c++
auto init = gasp2::initialize_catalysis(species, molar_masses,
                                        "input_gamma.xml");
if (!init) { /* handle error */ }
```

`initialize_catalysis` reads the XML file, checks the rules above, and caches the reaction set. The call expects species names and molar masses ordered consistently with later computations.

### Compute

```c++
auto result = gasp2::compute_catalysis_fluxes(T_wall, rho_wall);
if (!result) { /* handle error */ }
```

`compute_catalysis_fluxes` returns one surface flux per species in the input order for a wall temperature `T_wall` and near-wall mass densities `rho_wall`.

### Minimal Example

C++ program computing fluxes for five species:

```c++
#include <gasp2/gasp2.hpp>

int main() {
  std::vector<double> rho_wall{1.2, 1.4, 1.4, 0.1, 0.1};
  std::vector<std::string> species{"O2", "N2", "NO", "O", "N"};
  std::vector<double> molar_masses{32e-3, 28e-3, 30e-3, 16e-3, 14e-3};

  auto init = gasp2::initialize_catalysis(species, molar_masses,
                                          "input_gamma.xml");
  if (!init)
    return 1;

  auto fluxes = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!fluxes)
    return 1;
  for (double f : fluxes->cat_fluxes)
    std::cout << f << " ";
  std::cout << "\n";
  return 0;
}
```

Example XML (`inputs/input_gamma.xml`) with three GammaGiven reactions:

```xml
<gasp2 surface_model="catalysis">
  <reactions model="gamma_model">
    <reaction type="GammaGiven" formula="N + N -> N2">
      <gammas>N:0.5</gammas>
    </reaction>
    <reaction type="GammaGiven" formula="O + O -> O2">
      <gammas>O:0.5</gammas>
    </reaction>
    <reaction type="GammaGiven" formula="N + O -> NO">
      <gammas>N:0.5, O:0.5</gammas>
    </reaction>
  </reactions>
</gasp2>
```

## Finite-Rate Model

### Input File

Setting `<reactions model="finite_rate">` activates finite-rate surface kinetics. The XML file must include:

1. `<gasp2 surface_model="catalysis">` root element.
2. `<site_density>` specifying adsorption sites (m⁻²).
3. A series of `<reaction>` blocks, each declaring its `type` and reactant/product `formula`. Reaction types are case-insensitive (`Adsorption`, `Desorption`, `Adsorption/Desorption`, `Eley-Rideal` or `ER`, `Langmuir-Hinshelwood` or `LH`).

Parameters for each reaction are supplied inside a `<parameters>` element as
`name:value` pairs separated by commas or semicolons. Reversible
adsorption/desorption additionally uses a `<Kc>` block with the same
`name:value` format to provide equilibrium constants.

### Supported Reactions and Parameters

* **Adsorption** – gas species adsorb onto vacant sites `(s)`.
  * Parameters: `S0` (sticking coefficient at 0 K), `beta`, `Ead`.
* **Desorption** – adsorbed species desorb and free a site.
  * Parameters: `Ades` (pre-exponential factor), `beta`, `v` (vibrational
    frequency), `Edes`.
* **Adsorption/Desorption** – reversible adsorption plus desorption; include a
  `<Kc>` block specifying `A_eq`, `beta`, `K0`, and `dE`.
* **Eley-Rideal (ER)** – a gas species reacts directly with an adsorbed species.
  * Parameters: `gamma_er`, `beta`, `Eer`.
* **Langmuir-Hinshelwood (LH)** – two adsorbed species react on the surface.
  * Parameters: `C_lh`, `beta`, `Elh`.
Each reaction must be element-wise balanced. The library checks that all referenced species appear in the global species list and that the necessary parameters are present.

### Initialize

Initialization is identical to the gamma model but references a finite-rate input file:

```c++
auto init = gasp2::initialize_catalysis(species, molar_masses,
                                        "input_finite_rate.xml");
if (!init) { /* handle error */ }
```

### Compute

Flux evaluation uses the same API:

```c++
auto fluxes = gasp2::compute_catalysis_fluxes(T_wall, rho_wall);
if (!fluxes) { /* handle error */ }
```

### Minimal Example

C++ code:

```c++
#include <gasp2/gasp2.hpp>

int main() {
  std::vector<double> rho_wall{1.2, 1.4, 1.4, 0.1, 0.1};
  std::vector<std::string> species{"O2", "N2", "NO", "O", "N"};
  std::vector<double> molar_masses{32e-3, 28e-3, 30e-3, 16e-3, 14e-3};

  auto init = gasp2::initialize_catalysis(species, molar_masses,
                                          "input_finite_rate.xml");
  if (!init)
    return 1;

  auto fluxes = gasp2::compute_catalysis_fluxes(300.0, rho_wall);
  if (!fluxes)
    return 1;
  for (double f : fluxes->cat_fluxes)
    std::cout << f << " ";
  std::cout << "\n";
  return 0;
}
```

Sample XML highlighting adsorption and desorption of O using the expected
`<parameters>` format:

```xml
<gasp2 surface_model="catalysis">
  <reactions model="finite_rate">
    <site_density>1.0e-5</site_density>
    <reaction type="Adsorption" formula="O + (s) -> O(s)">
      <parameters>S0:0.3,beta:0.0,Ead:0.0</parameters>
    </reaction>
    <reaction type="Desorption" formula="O(s) -> O + (s)">
      <parameters>Ades:1.0e13,beta:0.0,v:1.0,Edes:4.0e5</parameters>
    </reaction>
  </reactions>
</gasp2>
```

An extended example demonstrating all reaction styles is provided at `inputs/finite_rate/input_typical_simple_finite_rate.xml`.
