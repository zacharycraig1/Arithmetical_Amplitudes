# Makefile for Arithmetic CHY Scattering verification
#
# Usage:
#   make verify          # Run all verification (requires Docker + wolframscript)
#   make verify-sage     # Run only Sage verification (requires Docker)
#   make verify-math     # Run only Mathematica verification (requires wolframscript)
#   make checksums       # Verify data integrity
#   make docker-pull     # Pull required Docker images

SAGE_CMD ?= sage
USE_DOCKER ?= 1
SAGE_IMAGE = sagemath/sagemath:10.4
# Cross-platform repo root for Docker bind mount.
# - On Unix, `pwd` is fine.
# - On Windows (Git Bash / MSYS / some make builds), `pwd` may not exist.
ifeq ($(OS),Windows_NT)
  DOCKER_MOUNT := $(CURDIR)
else
  DOCKER_MOUNT := $(shell pwd)
endif
DOCKER_RUN = docker run --rm -v "$(DOCKER_MOUNT):/home/sage/work" -w /home/sage/work $(SAGE_IMAGE)

# Detect wolframscript (optional)
ifeq ($(OS),Windows_NT)
  HAVE_WOLFRAM := $(shell where wolframscript >NUL 2>&1 && echo 1 || echo 0)
else
  HAVE_WOLFRAM := $(shell command -v wolframscript >/dev/null 2>&1 && echo 1 || echo 0)
endif

ifeq ($(USE_DOCKER),1)
  SAGE_RUN = $(DOCKER_RUN) sage
else
  SAGE_RUN = $(SAGE_CMD)
endif

# Portable directory creation for logs/
ifeq ($(OS),Windows_NT)
  MKDIR_LOGS = @if not exist logs mkdir logs
else
  MKDIR_LOGS = @mkdir -p logs
endif

.PHONY: all verify verify-full verify-big verify-sage verify-math verify-quick checksums docker-pull frob-quick clean test-8d11d chi8-mandelstam chi8-11d chi8-4d chi8-full frob-generic frob-11d verify-generic-quick

all: verify

# Pull Docker images
docker-pull:
	docker pull $(SAGE_IMAGE)

# Data integrity check
checksums:
	python tools/checksums.py data/SHA256SUMS

# Full verification
verify: verify-full
	@echo ""
	@echo "=================================="
	@echo "ALL VERIFICATIONS COMPLETE"
	@echo "=================================="

# Publication-grade verification (includes full n=7 run)
verify-full: verify-sage verify-n7-full
	@echo ""
ifeq ($(HAVE_WOLFRAM),1)
	$(MAKE) verify-math
else
	@echo "Skipping verify-math (wolframscript not found)."
endif
	@echo ""
	@echo "=================================="
	@echo "FULL VERIFICATION COMPLETE"
	@echo "=================================="

# Sage verification (Docker)
verify-sage: verify-n5-discriminant verify-n7-theorem

verify-n5-discriminant:
	@echo ""
	@echo "=== n=5 Discriminant Verification ==="
	$(SAGE_RUN) referee_checks/verify_n5_discriminant.sage --seeds 10

verify-n7-theorem:
	@echo ""
	@echo "=== n=7 Inert Prime Verification ==="
	$(MKDIR_LOGS)
	$(SAGE_RUN) referee_checks/verify_n7_theorem.sage --quick --method groebner --output logs/n7_quick.json

verify-n7-full:
	@echo ""
	@echo "=== n=7 Full Verification (slow) ==="
	$(MKDIR_LOGS)
	$(SAGE_RUN) referee_checks/verify_n7_theorem.sage --seeds 30 --method groebner --output logs/n7_full.json

verify-np2-extension:
	@echo ""
	@echo "=== n=7 Quadratic Extension Check ==="
	$(SAGE_RUN) referee_checks/verify_n7_np2_extension.sage --seed 0 --prime 7

verify-n7-np2: verify-np2-extension

verify-big: verify-n7-np2
	@echo ""
	@echo "=================================="
	@echo "BIG COMPUTE STEP COMPLETE"
	@echo "=================================="

# Mathematica verification
verify-math: verify-gram-levi

verify-gram-levi:
	@echo ""
	@echo "=== Gram-Levi-Civita Identity (Mathematica) ==="
	wolframscript -file referee_checks/verify_gram_levi.wl

# Quick verification (minimal runtime)
verify-quick:
	@echo ""
	@echo "=== Quick Verification ==="
	$(SAGE_RUN) referee_checks/verify_n5_discriminant.sage --seeds 3
	$(MKDIR_LOGS)
	$(SAGE_RUN) referee_checks/verify_n7_theorem.sage --quick --method groebner --output logs/n7_quick.json

# Quick Frobenius cycle-type extraction (small prime window)
frob-quick:
	@echo ""
	@echo "=== Frobenius Cycle Types (Quick) ==="
	$(MKDIR_LOGS)
	$(SAGE_RUN) code/sage/frob_cycletypes_eliminant.sage --seed_id 0 --seed_value 42 --p_min 5 --p_max 43 --outfile logs/frob_cycletypes_seed0_p5_43.json

# =============================================================================
# 8D/11D Generic Kinematics + χ₈ Correlation Tests
# =============================================================================

# Test the 8D/11D kinematics generator
test-8d11d:
	@echo ""
	@echo "=== 8D/11D Kinematics Smoke Test ==="
	cd code/sage && $(SAGE_RUN) test_8d11d_kinematics.sage

# χ₈ correlation test with dimension-free Mandelstams (fastest)
chi8-mandelstam:
	@echo ""
	@echo "=== χ₈ Correlation Test (Mandelstam mode, p≤100) ==="
	$(MKDIR_LOGS)
	cd code/sage && $(SAGE_RUN) chi8_correlation_test.sage --mode mandelstam --p_max 100 --outdir ../../logs

# χ₈ correlation test with 11D kinematics
chi8-11d:
	@echo ""
	@echo "=== χ₈ Correlation Test (11D mode, p≤100) ==="
	$(MKDIR_LOGS)
	cd code/sage && $(SAGE_RUN) chi8_correlation_test.sage --mode 8d11d --dim 11 --p_max 100 --outdir ../../logs

# χ₈ correlation test with 4D kinematics (comparison)
chi8-4d:
	@echo ""
	@echo "=== χ₈ Correlation Test (4D mode, p≤100) ==="
	$(MKDIR_LOGS)
	cd code/sage && $(SAGE_RUN) chi8_correlation_test.sage --mode 4d --p_max 100 --outdir ../../logs

# Full χ₈ sweep (all modes)
chi8-full:
	@echo ""
	@echo "=== Full χ₈ Correlation Sweep (all modes) ==="
	$(MKDIR_LOGS)
	cd code/sage && $(SAGE_RUN) chi8_correlation_test.sage --mode mandelstam --p_max 200 --outdir ../../logs
	cd code/sage && $(SAGE_RUN) chi8_correlation_test.sage --mode 8d11d --dim 11 --p_max 100 --outdir ../../logs

# Frobenius cycle types with generic kinematics
frob-generic:
	@echo ""
	@echo "=== Frobenius Cycle Types (Generic Kinematics) ==="
	$(MKDIR_LOGS)
	cd code/sage && $(SAGE_RUN) frob_cycletypes_generic.sage --kin_mode mandelstam --p_max 200 --outdir ../../logs

# Frobenius with 11D kinematics
frob-11d:
	@echo ""
	@echo "=== Frobenius Cycle Types (11D Kinematics) ==="
	$(MKDIR_LOGS)
	cd code/sage && $(SAGE_RUN) frob_cycletypes_generic.sage --kin_mode 8d11d --dim 11 --p_max 100 --outdir ../../logs

# Quick verification of new 8D/11D tools
verify-generic-quick:
	@echo ""
	@echo "=== Quick Verification of Generic Kinematics ==="
	cd code/sage && $(SAGE_RUN) test_8d11d_kinematics.sage
	$(MKDIR_LOGS)
	cd code/sage && $(SAGE_RUN) chi8_correlation_test.sage --mode mandelstam --p_max 50 --outdir ../../logs

# Clean generated files
clean:
	rm -f *.pyc *.sage.py
	rm -rf __pycache__
