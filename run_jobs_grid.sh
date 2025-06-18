#!/usr/bin/env bash

# grid_search.sh
# Usage: ./grid_search.sh GRID_TYPE
# GRID_TYPE choices:
#   glgt_theta  – loop over GLGT ∈ [0,0.4,0.8,2,10,100] and thetaLT ∈ [-3..3]
#   mu          – loop over Re(mu), Im(mu) ∈ [0,0.25,0.5,0.75,1.0]

if [ $# -lt 1 ]; then
  echo "Usage: $0 {glgt_theta|mu}"
  exit 1
fi

GRID_TYPE=$1
# common settings
N_EVENTS=100000
JOBS=1
DATE_LABEL="06_17_2025_100k"
BASE_OUTDIR="out"

case "$GRID_TYPE" in

  glgt_theta)
    glgt_values=(0 0.4 0.8 2 10 100)
    theta_values=(-3 -2 -1 0 1 2 3)
    for glgt in "${glgt_values[@]}"; do
      for theta in "${theta_values[@]}"; do
        OUTDIR="${DATE_LABEL}_GLGT=${glgt}_thetaLT=${theta}"
        echo "==== Running grid point: GLGT=${glgt}, thetaLT=${theta} ===="
        echo "Output directory: ${BASE_OUTDIR}/${OUTDIR}"

        ./run_jobs.rb \
          -n ${N_EVENTS} \
          -o "${OUTDIR}" \
          -j ${JOBS} \
          --no-gemc \
          --set "'StringSpinner:GLGT'=${glgt}" \
          --set "'StringSpinner:thetaLT'=${theta}" \
          --slurm

        echo
      done
    done
    ;;

  mu)
    # Note: the backslashes around (mu) are literal in the parameter name.
    re_mu_values=(0.1 0.5 0.8)
    im_mu_values=(0.1 0.5 0.8)
    for re_mu in "${re_mu_values[@]}"; do
      for im_mu in "${im_mu_values[@]}"; do
        OUTDIR="${DATE_LABEL}_ReMu=${re_mu}_ImMu=${im_mu}"
        echo "==== Running grid point: Re(μ)=${re_mu}, Im(μ)=${im_mu} ===="
        echo "Output directory: ${BASE_OUTDIR}/${OUTDIR}"

        ./run_jobs.rb \
          -n ${N_EVENTS} \
          -o "${OUTDIR}" \
          -j ${JOBS} \
          --no-gemc \
          --set "'StringSpinner:Re(mu)'=${re_mu}" \
          --set "'StringSpinner:Im(mu)'=${im_mu}" \


        echo
      done
    done
    ;;

  *)
    echo "Error: unknown GRID_TYPE '$GRID_TYPE'. Use 'glgt_theta' or 'mu'."
    exit 1
    ;;
esac
