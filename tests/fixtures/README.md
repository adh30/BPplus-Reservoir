# Test fixtures

These six XML files are real BP+/CardioScope device recordings, used as
regression-test fixtures for `read_BPplus.m`, `read_BPplusBPplus.m`, and
`read_BPplusCardioScope.m` in `../testReadBPplusReaders.m`.

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
