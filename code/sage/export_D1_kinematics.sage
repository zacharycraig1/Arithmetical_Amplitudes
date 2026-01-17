"""
Export D1 Kinematics to JSON
=============================

This script generates the explicit kinematic data for dataset D1
and exports it to JSON format for reproducibility.
"""

import json
import hashlib

load("kinematics.sage")

def export_D1_kinematics(output_file="../../data/D1_kinematics.json"):
    """Generate and export all D1 kinematics to JSON."""
    
    # Load seed specification
    with open("../../data/D1_seeds.json", "r") as f:
        seeds_data = json.load(f)
    
    results = {
        "dataset": "D1",
        "description": "Explicit 4D null momenta for n=7 CHY verification",
        "generation_method": "Pythagorean null vector parametrization with momentum conservation",
        "n_particles": 7,
        "kinematic_points": []
    }
    
    for seed_entry in seeds_data["seeds"]:
        seed = seed_entry["seed_value"]
        idx = seed_entry["id"]
        
        print(f"Generating seed {idx} (value={seed})...")
        
        try:
            momenta = make_4d_massless_point(7, seed)
            mandelstams = momenta_to_mandelstams(momenta)
            
            # Verification checks
            massless_check = all(check_null(p) for p in momenta)
            conservation_check = check_momentum_conservation(momenta)
            eps_value = levi_civita(momenta[0], momenta[1], momenta[2], momenta[3])
            genericity_check = eps_value != 0
            
            # Convert to serializable format
            momenta_list = []
            for a, p in enumerate(momenta):
                momenta_list.append({
                    "particle": a + 1,
                    "E": str(p[0]),
                    "x": str(p[1]),
                    "y": str(p[2]),
                    "z": str(p[3])
                })
            
            mandelstams_list = {}
            for (i, j), val in mandelstams.items():
                mandelstams_list[f"s_{i}{j}"] = str(val)
            
            point_data = {
                "id": idx,
                "seed": seed,
                "momenta": momenta_list,
                "mandelstams": mandelstams_list,
                "verification": {
                    "massless": massless_check,
                    "momentum_conservation": conservation_check,
                    "epsilon_nonzero": genericity_check,
                    "epsilon_value": str(eps_value)
                }
            }
            
            results["kinematic_points"].append(point_data)
            
        except Exception as e:
            print(f"  Failed: {e}")
            continue
    
    # Compute hash
    json_str = json.dumps(results, sort_keys=True, indent=2)
    sha256_hash = hashlib.sha256(json_str.encode()).hexdigest()
    results["sha256"] = sha256_hash
    
    # Write output
    with open(output_file, "w") as f:
        json.dump(results, f, indent=2)
    
    print(f"\nExported {len(results['kinematic_points'])} kinematic points")
    print(f"SHA256: {sha256_hash}")
    print(f"Output: {output_file}")
    
    return results


if __name__ == "__main__":
    export_D1_kinematics()
