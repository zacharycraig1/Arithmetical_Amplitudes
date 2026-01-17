"""
Query the transitive group database for 24T34.

Run (from publication_ready/):
  docker run --rm -v "$PWD:/home/sage/work" -w /home/sage/work sagemath/sagemath:10.4 sage publication_ready/tools/query_24T34.sage
"""

G = TransitiveGroup(24, 34)
print("TransitiveGroup(24,34)")
print("order =", G.order())
try:
    print("structure_description =", G.structure_description())
except Exception as e:
    print("structure_description failed:", e)
try:
    print("is_solvable =", G.is_solvable())
except Exception as e:
    print("is_solvable failed:", e)

# Also print GAP's StructureDescription directly (sometimes differs slightly)
try:
    print("gap StructureDescription =", gap(G).StructureDescription())
except Exception as e:
    print("gap StructureDescription failed:", e)

