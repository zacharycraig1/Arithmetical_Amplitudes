# Check transitivity/orbits for single-run tau+sigma pairs.
#
# Run:
#   docker run --rm -v "${PWD}:/home/sage/work" -w /home/sage/work sagemath/sagemath:10.4 gap -q tools/gap_check_single_runs.g

Print("=== Single-run checks ===\n\n");

# Shared tau (from logs)
tau := (1,2)(3,4)(5,6)(7,8)(9,10)(11,12)(13,14)(15,16)(17,18)(19,20)(21,22)(23,24);
Print("tau cycle structure: ", CycleStructurePerm(tau), "\n");
Print("tau orbits: ", Orbits(Group(tau), [1..24]), "\n\n");

# From D24_proof_run2.log
sigma_run2 := (3,7,5)(6,8);
G2 := Group(tau, sigma_run2);
Print("--- D24_proof_run2.log ---\n");
Print("sigma: ", sigma_run2, "\n");
Print("Order(G) = ", Size(G2), "\n");
Print("StructureDescription(G) = ", StructureDescription(G2), "\n");
Print("Transitive? ", IsTransitive(G2, [1..24]), "\n");
Print("Orbits: ", Orbits(G2, [1..24]), "\n\n");

# From D24_proof_run.log
sigma_run1 := (23,24);
G1 := Group(tau, sigma_run1);
Print("--- D24_proof_run.log ---\n");
Print("sigma: ", sigma_run1, "\n");
Print("Order(G) = ", Size(G1), "\n");
Print("StructureDescription(G) = ", StructureDescription(G1), "\n");
Print("Transitive? ", IsTransitive(G1, [1..24]), "\n");
Print("Orbits: ", Orbits(G1, [1..24]), "\n\n");

# From D24_disc_only_quick2.log
sigma_disc := (12,14);
Gd := Group(tau, sigma_disc);
Print("--- D24_disc_only_quick2.log ---\n");
Print("sigma: ", sigma_disc, "\n");
Print("Order(G) = ", Size(Gd), "\n");
Print("StructureDescription(G) = ", StructureDescription(Gd), "\n");
Print("Transitive? ", IsTransitive(Gd, [1..24]), "\n");
Print("Orbits: ", Orbits(Gd, [1..24]), "\n\n");

Print("=== Reference: 24T34 ===\n");
T := TransitiveGroup(24,34);
Print("Order(24T34) = ", Size(T), "\n");
Print("StructureDescription(24T34) = ", StructureDescription(T), "\n");
Print("Transitive? ", IsTransitive(T, [1..24]), "\n");
Print("Orbits: ", Orbits(T, [1..24]), "\n");

