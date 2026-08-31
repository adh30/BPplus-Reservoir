# Test fixtures

Two kinds of fixture live here.

The six `CardioScope_*.xml` / `BPplus_*.xml` files are real BP+/CardioScope
device recordings, used as regression-test fixtures for `read_BPplus.m`,
`read_BPplusBPplus.m`, and `read_BPplusCardioScope.m` in
`../testReadBPplusReaders.m`.

The three `synthetic_ai_v2_Type{A,B,C}.xml` files are NOT real recordings -
see "Synthetic ai_v2.m fixtures" below.

## Anonymisation status

- `datetime` on every file has been replaced with `2020-01-01T00:00:00`.
- No patient name, date of birth, address, or MRN fields exist in this
  XML schema, so there was nothing of that kind to remove.
- **`device_id`, `guid`, and `nibp_id` have NOT been anonymised.** These
  are hardware serial numbers / per-recording identifiers, not patient
  identifiers, but combined with a real-world appointment or
  device-assignment log they could act as a quasi-identifier for a
  specific visit. Worth reviewing before treating this repo as fully
  free of identifying information, particularly if it's ever made public
  and that matters for your use case.

## Naming

Files are named `<Type>_<version>[_nnn].xml`, where `<Type>` is the XML
root element (`BPplus` or `CardioScope`) and `<version>` is the
`firmware_version` attribute where present, else `software_version`. A
`_nnn` suffix disambiguates files that would otherwise share a name.

| File | Schema | Notes |
|---|---|---|
| `CardioScope_SW.R7.VME.035_001.xml` | CardioScope | older recording (2012) |
| `CardioScope_SW.R7.VME.035_002.xml` | CardioScope | newer recording (2013), same device |
| `BPplus_4.4.44094.43303_001.xml` | BPplus | has `cST`, `cAveragePulsePointsIndexes` |
| `BPplus_4.4.44094.43303_002.xml` | BPplus | has `cST`, `cAveragePulsePointsIndexes` |
| `BPplus_4.4.44094.43303_003.xml` | BPplus | has `cST`, `cAveragePulsePointsIndexes` |
| `BPplus_DeviceFirmwareVersion.xml` | BPplus | older/reduced schema; placeholder `device_id`/`firmware_version`; no `cST` or `cAveragePulsePointsIndexes` |

## Synthetic ai_v2.m fixtures

`synthetic_ai_v2_TypeA.xml`, `synthetic_ai_v2_TypeB.xml`, and
`synthetic_ai_v2_TypeC.xml` are used by `../testAiV2WaveformTypes.m` to
check `ai_v2.m`'s inflection-point detection and Murgo et al. (Circulation
1980;62(1):105-16) Type A/B/C classification. They are each a copy of the
real, already-anonymised `BPplus_4.4.44094.43303_001.xml`, with only the
following fields changed:

- `cAveragePulse` - replaced with a hand-constructed waveform (sum of two
  skewed-Gaussian components on a baseline, sampled at 200 Hz - the real
  device rate - over ~600 ms, matching a real fixture's average-beat
  length). This is the only field that reaches `ai_v2.m` (it becomes
  `ao.p_av`).
- `cSys` / `cDia` / `cPP` - updated to match the synthetic waveform's own
  SBP/DBP/PP, for cosmetic consistency only. Not used by `ai_v2.m`.
- `device_id`, `guid`, `nibp_id` - replaced with obvious placeholders
  (`SYNTHETIC0000Type{A,B,C}`, etc.), since this is not a real device
  recording.

Every other field (`sBaseLined`, `baEstimate`, `cEstimate`, pulse-start
indexes, the raw base64 waveform blobs, all the derived cAP/cAIx/cSEVR/etc.
clinical fields, ...) is an **untouched, stale copy** of the original real
recording. `read_BPplus.m` needs internally-consistent index/bounds data
in those fields to run without erroring, but none of them feed `ai_v2.m`,
and none of them describe the synthetic waveform - do not read any
clinical meaning into them for these three files.

Provenance of the waveform shapes and target classification values:

- Sample rate (200 Hz) and pulse duration (~600 ms): read from real
  fixtures, not assumed.
- Type A target (ai = +29, i.e. deltaP/PP = +0.29) and Type C target
  (ai = -19, deltaP/PP = -0.19): taken directly from Murgo JP, Westerhof N,
  Giolma JP, Altobelli SA. Aortic input impedance in normal man:
  relationship to pressure wave forms. Circulation 1980;62(1):105-16 -
  Figure 3's own worked numeric examples for a Type A and a Type C beat.
  Figure 3 itself carries no printed time-axis calibration and does not
  correspond exactly to any individual patient row in the paper's Table 1
  (checked: no patient has PP=57/deltaP=17 or PP=31/deltaP=-6), so it was
  treated as an illustrative schematic, not a traceable single-patient
  recording - only its shape and the two numeric ratios were used.
- Type B target (ai = +9): NOT sourced from any reference. Murgo 1980 and
  the related 1981 exercise paper (Murgo et al., Circ Res 1981;48:334-43)
  both describe Type B only as "intermediate between types A and C...not
  illustrated," bounded 0.0 < deltaP/PP < 0.12. +9 was chosen by this
  assistant to sit mid-band; treat it as a design choice, not a citable
  value.
- The two-skewed-Gaussian-plus-baseline waveform construction itself is a
  standard pulse-decomposition-style synthesis technique, not a formula
  given in Murgo, Kelly (Kelly R, Hayward C, Avolio A, O'Rourke M.
  Noninvasive determination of age-related changes in the human arterial
  pulse. Circulation 1989;80(6):1652-9), or Karamanoglu (Karamanoglu M. A
  system for analysis of arterial blood pressure waveforms in humans.
  Comput Biomed Res 1997;30:244-55). Those references informed the
  qualitative shape (shoulder/dip morphology, foot-to-peak structure) and
  the fiducial-point definitions ai_v2.m implements, not the parametric
  form used to generate these fixtures.
- Expected `typetxt`/`ai` values in `testAiV2WaveformTypes.m` come from a
  Python reimplementation of `ai_v2.m`, not from running MATLAB - see the
  caveat at the top of that test file.
