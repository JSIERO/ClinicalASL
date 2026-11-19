# ClinicalASL nightly verification (unattended, ~2 h). Writes G:\DATA\MOYAMOYA\IMAGER\tests\<Name>\summary.txt
# (no patient data: counts, exit codes, comparison tables, redacted dciodvfy summaries).
#   export : step 12 only on two copies of IMAGER\workingdir: the image's own code (reference) vs the code under test,
#            then compare_dicom (headers + pixels, UIDs/times ignored) and compare_maps
#   cvr    : full CVR run on ASLtest            single : full single-mode run on the pre-ACZ part of ASLtest
#   dciodvfy on both outputs (logs in IMAGER\dciodvfy\<Name>_*), compare_maps/compare_all/compare_dicom vs -Ref
# Code under test: the repo's working tree mounted at /fix (default) or the image's own code (-ImageCode).
# Schedule (one-off, tonight 01:00; PC awake, Docker Desktop running, user logged on):
#   $a = New-ScheduledTaskAction -Execute powershell.exe -Argument '-NoProfile -ExecutionPolicy Bypass -File "<repo>\tools\nightly_verify.ps1" -Name <name>'
#   Register-ScheduledTask ClinicalASL_nightly_<name> -Action $a -Trigger (New-ScheduledTaskTrigger -Once -At 01:00) -Settings (New-ScheduledTaskSettingsSet -WakeToRun -StartWhenAvailable)
param(
  [Parameter(Mandatory = $true)][string]$Name,
  [string]$Image = "jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager:1.2.0",
  [string]$Ref = "G:\DATA\MOYAMOYA\IMAGER\tests\v1.2.0_image",
  [string[]]$Steps = @("export", "cvr", "single", "dciodvfy", "compare"),
  [switch]$ImageCode
)
$ErrorActionPreference = "Continue"
$repo = Split-Path -Parent $PSScriptRoot
$tools = "$repo\tools"
$imager = "G:\DATA\MOYAMOYA\IMAGER"
$t = "$imager\tests\$Name"
New-Item -ItemType Directory -Force $t | Out-Null
$sum = "$t\summary.txt"
function Log($m) { "$(Get-Date -Format 'yyyy-MM-dd HH:mm') $m" | Out-File -Append -Encoding utf8 $sum }
function Run-Tool($title, $script) { Log "== $title"; & python "$tools\$script" @args 2>&1 | Out-File -Append -Encoding utf8 $sum }
$fix = if ($ImageCode) { @() } else { @("-w", "/fix", "-e", "PYTHONPATH=/fix", "-v", "$repo\python:/fix:ro") }
$py = if ($ImageCode) { "/opt/ClinicalASL/python" } else { "/fix" }
Log "START $Name image=$Image code=$(if ($ImageCode) {'image'} else {'working tree ' + (git -C $repo rev-parse --short HEAD) + ' + local changes'}) steps=$($Steps -join ',')"

function Check-Output($dir) {
  $dcm = @(Get-ChildItem "$dir\*.dcm").Count
  $err = @(Select-String "$dir\clinicalasl.log" -Pattern ' - ERROR - ' -ErrorAction SilentlyContinue).Count
  $pid_ = @(Select-String "$dir\clinicalasl.log" -Pattern 'PatientID' -ErrorAction SilentlyContinue).Count
  $lines = @(Get-Content "$dir\clinicalasl.log" -ErrorAction SilentlyContinue).Count
  Log "   $dir : DICOM files $dcm, log lines $lines, ERROR lines $err, PatientID lines $pid_"
}

function Run-Docker($label, $inDir, $out, $work, $script, $runArgs) {
  New-Item -ItemType Directory -Force $out, $work | Out-Null
  $mounts = @("-v", "${out}:/output", "-v", "${work}:/workingdir")
  if ($inDir) { $mounts += @("-v", "${inDir}:/input:ro") }
  $t0 = Get-Date
  & docker run --rm --network none --name "clinicalasl_$($Name)_$label" @fix @mounts --entrypoint python $Image "$py/$script" @runArgs *> "$t\$label`_docker_stdout.log"
  $code = $LASTEXITCODE
  Log "$label exit=$code minutes=$([math]::Round(((Get-Date) - $t0).TotalMinutes, 1))"
  Check-Output $out
}

if ($Steps -contains "export") {
  foreach ($v in @("export_ref", "export_new")) {
    robocopy "$imager\workingdir" "$t\$v\work" /E /NFL /NDL /NJH /NJS /NP | Out-Null
  }
  $saveFix = $fix; $savePy = $py
  $fix = @(); $py = "/opt/ClinicalASL/python"
  Run-Docker "export_ref" $null "$t\export_ref\output" "$t\export_ref\work" "export_results.py" @("/workingdir", "/output")
  $fix = $saveFix; $py = $savePy
  Run-Docker "export_new" $null "$t\export_new\output" "$t\export_new\work" "export_results.py" @("/workingdir", "/output")
  Run-Tool "export: DICOM reference code vs code under test" "compare_dicom.py" "$t\export_ref\output" "$t\export_new\output"
  Run-Tool "export: maps reference vs new (expect identical)" "compare_maps.py" "$t\export_ref\work\ASL" "$t\export_new\work\ASL"
  Remove-Item -Recurse -Force "$t\export_ref\work", "$t\export_new\work" -Confirm:$false  # patient DICOM copies
}
if ($Steps -contains "cvr") {
  Run-Docker "cvr" "$imager\ASLtest" "$t\cvr\output" "$t\cvr\work" "run_pipeline.py" @("/input", "/output", "/workingdir", "--mode", "cvr")
}
if ($Steps -contains "single") {
  & python "$tools\make_preacz_input.py" "$imager\ASLtest" "$t\single\input" | Out-Null
  Run-Docker "single" "$t\single\input" "$t\single\output" "$t\single\work" "run_pipeline.py" @("/input", "/output", "/workingdir", "--mode", "single")
  Remove-Item -Recurse -Force "$t\single\input" -Confirm:$false  # patient DICOM copy
}
foreach ($m in @("cvr", "single")) {
  if (-not (Test-Path "$t\$m\output")) { continue }
  # keep the evidence lean: no copies of source DICOMs / NIfTIs in the archived work folder
  Remove-Item -Recurse -Force "$t\$m\work\DICOMORIG", "$t\$m\work\NIFTI" -Confirm:$false -ErrorAction SilentlyContinue
  if ($Steps -contains "dciodvfy") { Run-Tool "dciodvfy $m" "run_dciodvfy.py" "$t\$m\output" "$imager\dciodvfy\$($Name)_$m" }
  if ($Steps -contains "compare" -and (Test-Path "$Ref\$m")) {
    Run-Tool "compare_maps $m vs $Ref" "compare_maps.py" "$Ref\$m\work\ASL" "$t\$m\work\ASL"
    Run-Tool "compare_dicom $m vs $Ref (pixels differ by ANTs noise; headers should match)" "compare_dicom.py" "$Ref\$m\output" "$t\$m\output"
  }
}
Log "ALL DONE"
