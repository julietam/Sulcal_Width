#!/bin/bash

# =========================================================================
# qsub Submission Script for Pipeline Stage B (Template)
# SGE Version with Dependency - ADAPT FOR YOUR SCHEDULER AND ENVIRONMENT
# =========================================================================
# Runs StageB_CalculateMetrics.m, assuming Stage A outputs exist.
# Waits for Stage A job array to complete successfully before starting.

# --- Job Submission Options (ADAPT FOR YOUR SGE SYSTEM or REPLACE FOR OTHER SCHEDULERS!) ---
#$ -N <Your_Job_Name_B>     # <<< Choose a Job Name for Stage B >>>
#$ -l h_rt=<Your_Time_Limit>       # <<< Request runtime PER TASK (e.g., 01:00:00 for 1 hour) - Often less than Stage A >>>
#$ -l h_vmem=<Your_Memory_Limit>   # <<< Request RAM PER TASK (e.g., 4G for 4 GB) - Often less than Stage A >>>
#$ -pe smp 1              # <<< Request 1 core per task (Syntax/option might vary) >>>
#$ -t 1-<N>                 # <<< SET N to the EXACT number of subjects in your subject_list.txt (MUST MATCH STAGE A) >>>

# >>> IMPORTANT: JOB DEPENDENCY <<<
# Replace <StageA_JobID> with the ACTUAL Job ID returned when you submitted Stage A
# Example: If Stage A job ID was 12345, use: #$ -hold_jid 12345
#$ -hold_jid <StageA_JobID>  # <<< EDIT THIS MANUALLY AFTER SUBMITTING STAGE A! Syntax is SGE specific >>>

#$ -o /path/to/your/logs/pipeline_B_job_$TASK_ID.out # <<< Set path for stdout logs. Ensure directory exists and is writable! >>>
#$ -e /path/to/your/logs/pipeline_B_job_$TASK_ID.err # <<< Set path for stderr logs. Ensure directory exists and is writable! >>>
#$ -cwd                   # Run job from current directory (Optional but often useful)
#$ -m n                   # No email notifications (or 'e' for end, 'a' for abort, 'b' for begin)
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
# !!! IMPORTANT: Set these paths correctly (should match Stage A script) !!!
SCRIPT_DIR="/path/to/your/scripts_and_subject_list" # <<< Directory containing .m scripts AND subject_list.txt >>>
SUBJECT_LIST="${SCRIPT_DIR}/subject_list.txt"    # Assumes subject_list.txt is in SCRIPT_DIR
MATLAB_SCRIPT_B="StageB_CalculateMetrics"      # Name of Stage B script (without .m)
LOG_DIR="/path/to/your/logs"                   # <<< Directory for MATLAB specific logs (matlab_B_*.log) >>>

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

echo "===== Processing Subject (Stage B): ${SUBJECT_ID} (Task ID: ${TASK_ID}) ====="

# --- Export Subject ID for MATLAB ---
# StageB_CalculateMetrics.m reads this environment variable
export JOB_SUBJECT_ID=${SUBJECT_ID}

# --- Define MATLAB options ---
# -batch runs script non-interactively and exits; -nodisplay/-nosplash prevent GUI elements
MATLAB_OPTS="-nodisplay -nosplash -batch"

# --- Run Stage B ---
echo "$(date): Running Stage B for ${SUBJECT_ID}..."
MATLAB_LOG_B="${LOG_DIR}/matlab_B_${SUBJECT_ID}.log"
# Construct MATLAB command: Change directory, run script
# Ensure StageB_CalculateMetrics.m exists in SCRIPT_DIR
CMD_B="cd('${SCRIPT_DIR}'); ${MATLAB_SCRIPT_B}" # Simpler call assuming the script handles its own try/catch/exit

echo "Executing MATLAB Command: ${MATLAB_CMD} ${MATLAB_OPTS} \"${CMD_B}\""

# Execute MATLAB - redirect stdout and stderr from MATLAB execution to the log file
${MATLAB_CMD} ${MATLAB_OPTS} "${CMD_B}" > ${MATLAB_LOG_B} 2>&1
STATUS_B=$? # Capture MATLAB's exit status

# Check Stage B exit status
if [ $STATUS_B -ne 0 ]; then
  # Send error message to the job's standard error log
  echo "$(date): Error: MATLAB Stage B failed for ${SUBJECT_ID}. Status: ${STATUS_B}. Check ${MATLAB_LOG_B} and job error log (${LOG_DIR}/pipeline_B_job_${TASK_ID}.err)." >&2
  exit $STATUS_B # Exit the job script with failure status
fi

echo "$(date): Stage B completed successfully for ${SUBJECT_ID}."
echo "-----------------------------------------------------"
exit 0
# --- End of Script ---
