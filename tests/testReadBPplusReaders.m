classdef testReadBPplusReaders < matlab.unittest.TestCase
    % Regression tests for read_BPplus.m / read_BPplusBPplus.m /
    % read_BPplusCardioScope.m against real recordings.
    %
    % Fixtures in tests/fixtures/ are real device output files with the
    % datetime attribute replaced (2020-01-01T00:00:00) and renamed by
    % schema type + firmware version. See tests/fixtures/README.md for
    % provenance and residual identifiers (device_id/guid/nibp_id are NOT
    % anonymised).
    %
    % These tests encode values that were manually verified against the
    % raw XML text when the beta7 fixes were made (see version history in
    % bpp_Res2.m). They have NOT been executed against MATLAB - review
    % and run before relying on them, and report back anything that
    % fails so the assertions (or the code) can be corrected.

    properties (Constant)
        Npoly = 3;
        Frame = 9;
    end

    properties
        FixtureDir
    end

    methods (TestClassSetup)
        function addRepoToPath(testCase)
            here = fileparts(mfilename('fullpath'));
            repoRoot = fileparts(here);
            addpath(repoRoot);
            testCase.FixtureDir = [fullfile(here, 'fixtures') filesep];
        end
    end

    methods (Test)

        function testCardioScope001_DispatchAndBasics(testCase)
            % originally MeasData_00017.xml (2012 recording)
            [data, metadata, ~, ba, ao] = read_BPplus(testCase.FixtureDir, ...
                'CardioScope_SW.R7.VME.035_001.xml', testCase.Npoly, testCase.Frame);

            testCase.verifyTrue(isfield(data, 'CardioScope'));
            testCase.verifyEqual(ba.sbp, 106);
            testCase.verifyEqual(ba.dbp, 68);
            testCase.verifyEqual(metadata.samplerate, 200);
            testCase.verifyEqual(metadata.snr, 21);

            % software_version "SW.R7.VME.035" -> parsed build 35 -> not
            % build 38 -> default lag
            testCase.verifyEqual(ao.lag, 0.18);

            % CardioScope schema has no central duration-of-systole field
            testCase.verifyTrue(isnan(ao.ed));

            testCase.verifyTrue(isfield(metadata, 'patient_id'));
            testCase.verifyFalse(isfield(metadata, 'patienet_id'));
        end

        function testCardioScope002_DispatchAndBasics(testCase)
            % originally Cardioscope_00480.xml (2013 recording)
            [data, metadata, ~, ba, ao] = read_BPplus(testCase.FixtureDir, ...
                'CardioScope_SW.R7.VME.035_002.xml', testCase.Npoly, testCase.Frame);

            testCase.verifyTrue(isfield(data, 'CardioScope'));
            testCase.verifyEqual(ba.sbp, 123);
            testCase.verifyEqual(ba.dbp, 79);
            testCase.verifyEqual(metadata.samplerate, 200);
            testCase.verifyEqual(metadata.snr, 23);
            testCase.verifyEqual(ao.lag, 0.18);
            testCase.verifyTrue(isnan(ao.ed));
        end

        function testBPplus001_AveragePulsePointsIndexesRegression(testCase)
            % originally new_test3.xml
            [data, metadata, ~, ba, ao] = read_BPplus(testCase.FixtureDir, ...
                'BPplus_4.4.44094.43303_001.xml', testCase.Npoly, testCase.Frame);

            testCase.verifyTrue(isfield(data, 'BPplus'));
            testCase.verifyEqual(ba.sbp, 111);
            testCase.verifyEqual(ba.dbp, 69);
            testCase.verifyEqual(metadata.snr, 10);

            % Regression test for the duplicate-block bug (beta7): must
            % come from cAveragePulsePointsIndexes ("0,10,28,40,50,117",
            % 6 values) not cPulseStartIndexes (17 values). The old buggy
            % code would have produced a 17-element array here instead.
            testCase.verifyEqual(ao.averagePulsePointsIndexes, [1;11;29;41;51;118]);
            testCase.verifyEqual(numel(ao.averagePulsePointsIndexes), 6);

            % Regression test for the NIBP mode guard fix (beta7): this
            % field never got set under the old buggy isfield check.
            testCase.verifyTrue(isfield(metadata, 'mode'));
            testCase.verifyEqual(metadata.mode, "MeasureOnDeflate");

            % ao.ed from cST=255ms -> 0.255s
            testCase.verifyEqual(ao.ed, 0.255, 'AbsTol', 1e-9);

            testCase.verifyTrue(isfield(metadata, 'patient_id'));
            testCase.verifyFalse(isfield(metadata, 'patienet_id'));
        end

        function testBPplus002_AveragePulsePointsIndexesRegression(testCase)
            % originally new_test1.xml
            [data, metadata, ~, ba, ao] = read_BPplus(testCase.FixtureDir, ...
                'BPplus_4.4.44094.43303_002.xml', testCase.Npoly, testCase.Frame);

            testCase.verifyTrue(isfield(data, 'BPplus'));
            testCase.verifyEqual(ba.sbp, 128);
            testCase.verifyEqual(ba.dbp, 75);
            testCase.verifyEqual(metadata.snr, 11);

            % cAveragePulsePointsIndexes = "0,10,29,46,63,279" (6 values)
            testCase.verifyEqual(ao.averagePulsePointsIndexes, [1;11;30;47;64;280]);
            testCase.verifyEqual(numel(ao.averagePulsePointsIndexes), 6);

            testCase.verifyTrue(isfield(metadata, 'mode'));
            testCase.verifyEqual(metadata.mode, "MeasureOnDeflate");

            % ao.ed from cST=320ms -> 0.320s
            testCase.verifyEqual(ao.ed, 0.320, 'AbsTol', 1e-9);
        end

        function testBPplus003_AveragePulsePointsIndexesRegression(testCase)
            % originally new_test2.xml
            [data, metadata, ~, ba, ao] = read_BPplus(testCase.FixtureDir, ...
                'BPplus_4.4.44094.43303_003.xml', testCase.Npoly, testCase.Frame);

            testCase.verifyTrue(isfield(data, 'BPplus'));
            testCase.verifyEqual(ba.sbp, 125);
            testCase.verifyEqual(ba.dbp, 70);
            testCase.verifyEqual(metadata.snr, 13);

            % cAveragePulsePointsIndexes = "0,10,28,44,80,172" (6 values)
            testCase.verifyEqual(ao.averagePulsePointsIndexes, [1;11;29;45;81;173]);
            testCase.verifyEqual(numel(ao.averagePulsePointsIndexes), 6);

            testCase.verifyTrue(isfield(metadata, 'mode'));
            testCase.verifyEqual(metadata.mode, "MeasureOnDeflate");

            % ao.ed from cST=405ms -> 0.405s
            testCase.verifyEqual(ao.ed, 0.405, 'AbsTol', 1e-9);
        end

        function testBPplusDeviceFirmwareVersion_MissingFieldsHandledSafely(testCase)
            % originally old_test1.xml - older/reduced schema: no cST, no
            % cAveragePulsePointsIndexes at all
            [data, metadata, ~, ba, ao] = read_BPplus(testCase.FixtureDir, ...
                'BPplus_DeviceFirmwareVersion.xml', testCase.Npoly, testCase.Frame);

            testCase.verifyTrue(isfield(data, 'BPplus'));
            testCase.verifyEqual(ba.sbp, 168);
            testCase.verifyEqual(ba.dbp, 88);
            testCase.verifyEqual(metadata.snr, 22);

            % cAveragePulsePointsIndexes absent from this file's schema -
            % the isfield guard means the field should simply never be set
            testCase.verifyFalse(isfield(ao, 'averagePulsePointsIndexes'));

            % cST absent -> ao.ed must be NaN, not -1 or []
            testCase.verifyTrue(isnan(ao.ed));

            % NibpModeUsed = "Undefined" is still present as a field, so
            % metadata.mode should still be set (to the literal text
            % "Undefined", which is the device's own value, not a MATLAB
            % missing-value marker)
            testCase.verifyTrue(isfield(metadata, 'mode'));
            testCase.verifyEqual(metadata.mode, "Undefined");
        end

    end
end
