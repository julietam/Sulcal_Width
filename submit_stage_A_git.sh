#!/bin/bash

# =========================================================================
# qsub Submission Script for Pipeline Stage A (Template)
# SGE Version - ADAPT FOR YOUR SCHEDULER AND ENVIRONMENT
# =========================================================================

# --- Job Submission Options (ADAPT FOR YOUR SGE SYSTEM or REPLACE FOR OTHER SCHEDULERS!) ---
#$ -N <Your_Job_Name_A>     # <<< Choose a Job Name for Stage A >>>
#$ -l h_rt=<Your_Time_Limit>       # <<< Request runtime PER TASK (e.g., 02:00:00 for 2 hours) >>>
#$ -l h_vmem=<Your_Memory_Limit>   # <<< Request RAM PER TASK (e.g., 8G for 8 GB) >>>
#$ -pe smp 1             # <<< Request 1 core per task (Syntax/option might vary) >>>
#$ -t 1-<N>                # <<< SET N to the EXACT number of subjects in your subject_list.txt >>>
#$ -o /path/to/your/logs/pipeline_A_job_$TASK_ID.out # <<< Set path for stdout logs. Ensure directory exists and is writable! >>>
#$ -e /path/to/your/logs/pipeline_A_job_$TASK_ID.err # <<< Set path for stderr logs. Ensure directory exists and is writable! >>>
#$ -cwd                  # Run job from current directory (Optional but often useful)
#$ -m n                  # No email notifications (or 'e' for end, 'a' for abort, 'b' for begin)
#$ -M <your_email@example.com> # <<< Add your email if using notifications (-m option) >>>

# --- Environment Setup ---
echo "-----------------------------------------------------"
echo "Job Started: $(date)"
echo "Running on host: $(hostname)"
# This line assumes SGE which provides $SGE_TASK_ID. Other schedulers use different variables (e.g., $SLURM_ARRAY_TASK_ID)
echo "Scheduler Task ID: $SGE_TASK_ID" # <<< Adapt if not using SGE >>>

# Load MATLAB module if necessary (Uncomment and adapt for your cluster)
# module load matlab/<Your_Matlab_Version> # <<< Specify your MATLAB version/module name >>>

# --- Define Paths ---
# !!! IMPORTANT: Set these paths correctly !!!
SCRIPT_DIR="/path/to/your/scripts_and_subject_list" # <<< Directory containing .m scripts AND subject_list.txt >>>
SUBJECT_LIST="${SCRIPT_DIR}/subject_list.txt"    # Assumes subject_list.txt is in SCRIPT_DIR
MATLAB_SCRIPT_A="StageA_GenerateVolumes2"        # Name of Stage A script (without .m)
LOG_DIR="/path/to/your/logs"                   # <<< Directory for MATLAB specific logs (matlab_A_*.log) >>>

# !!! IMPORTANT: Set MATLAB_CMD to the full path from 'which matlab' on compute nodes !!!
MATLAB_CMD="/path/to/your/matlab"             # <<< Full path to MATLAB executable >>>

# --- Create Log Directory ---
# Note: The directories for SGE logs (-o, -e) and MATLAB logs (LOG_DIR) should exist beforehand
mkdir -p ${LOG_DIR}
if [ ! -d "${LOG_DIR}" ] || [ ! -w "${LOG_DIR}" ]; then
    echo "Error: MATLAB log directory ${LOG_DIR} does not exist or is not writable." >&2
    exit 1
fi
echo "MATLAB Log directory: ${LOG_DIR}"

# --- Get Subject ID for this Task ---
# Assumes SGE ($SGE_TASK_ID). Adapt variable if using a different scheduler.
TASK_ID=$SGE_TASK_ID # <<< Adapt if not using SGE >>>
if [ -z "$TASK_ID" ]; then
    echo "Error: Scheduler Task ID environment variable is not set." >&2
    exit 1
fi

# Extract the Nth subject ID from the list file (N = TASK_ID)
echo "Subject list file: ${SUBJECT_LIST}"
if [ ! -f "${SUBJECT_LIST}" ]; then
    echo "Error: Subject list file not found: ${SUBJECT_LIST}" >&2
    exit 1
fi
SUBJECT_ID=$(sed -n "${TASK_ID}p" ${SUBJECT_LIST})
if [ -z "$SUBJECT_ID" ]; then
    echo "Error: No Subject ID found for Task ${TASK_ID} in ${SUBJECT_LIST}" >&2
    # Exit cleanly if task ID > number of subjects listed (optional, depends on scheduler behavior)
    echo "Exiting successfully as task ID might be beyond subject list size."
    exit 0
fi

echo "Processing Subject (Stage A): ${SUBJECT_ID} (Task ID: ${TASK_ID})"

# --- Export Subject ID for MATLAB ---
# StageA_GenerateVolumes2.m reads this environment variable
export JOB_SUBJECT_ID=${SUBJECT_ID}

# --- Define MATLAB options ---
# -batch runs script non-interactively and exits; -nodisplay/-nosplash prevent GUI elements
MATLAB_OPTS="-nodisplay -nosplash -batch"

# --- Run Stage A ---
echo "$(date): Running Stage A for ${SUBJECT_ID}..."
echo "Script Directory: ${SCRIPT_DIR}"
echo "MATLAB Command: ${MATLAB_CMD}"
MATLAB_LOG_A="${LOG_DIR}/matlab_A_${SUBJECT_ID}.log"

# Construct MATLAB command: Change directory, run script, handle errors, exit
# Ensures script runs from the correct directory
CMD_A="cd('${SCRIPT_DIR}'); ${MATLAB_SCRIPT_A}" # Simpler call assuming the script handles its own try/catch/exit

echo "Executing MATLAB Command: ${MATLAB_CMD} ${MATLAB_OPTS} \"${CMD_A}\""

# Execute MATLAB - redirect stdout and stderr from MATLAB execution to the log file
${MATLAB_CMD} ${MATLAB_OPTS} "${CMD_A}" > ${MATLAB_LOG_A} 2>&1
STATUS_A=$? # Capture MATLAB's exit status

# Check Stage A exit status
if [ $STATUS_A -ne 0 ]; then
  # Send error message to the job's standard error log
  echo "$(date): Error: MATLAB Stage A failed for ${SUBJECT_ID}. Status: ${STATUS_A}. Check ${MATLAB_LOG_A} and job error log (${LOG_DIR}/pipeline_A_job_${TASK_ID}.err)." >&2
  exit $STATUS_A # Exit the job script with failure status
fi

echo "$(date): Stage A completed successfully for ${SUBJECT_ID}."
echo "-----------------------------------------------------"
exit 0
# --- End of Script ---
