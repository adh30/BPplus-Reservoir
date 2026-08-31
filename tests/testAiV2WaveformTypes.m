classdef testAiV2WaveformTypes < matlab.unittest.TestCase
    % Regression/validation tests for ai_v2.m's inflection-point detection
    % and Murgo et al. (Circulation 1980;62(1):105-16) Type A/B/C
    % classification, using synthetic aortic average-beat waveforms rather
    % than a real recording.
    %
    % Fixtures in tests/fixtures/synthetic_ai_v2_Type{A,B,C}.xml are NOT
    % real device recordings. Each is a real, already-anonymised BPplus
    % fixture (BPplus_4.4.44094.43303_001.xml) with only cAveragePulse
    % (and, cosmetically, cSys/cDia/cPP) replaced by a hand-constructed
    % waveform designed to fall into that Murgo type; device_id/guid/
    % nibp_id were also replaced with obvious placeholders. Every other
    % field (sBaseLined, baEstimate, cEstimate, pulse-start indexes, raw
    % waveform blobs, etc.) is an untouched copy of the real recording,
    % kept only so read_BPplus.m's unrelated index/bounds logic has
    % internally-consistent data to run against - none of it feeds
    % ai_v2.m. See tests/fixtures/README.md for the full provenance and
    % the assumptions behind the synthetic waveform (sample rate and
    % pulse duration are taken from real fixtures; the two-component
    % pulse shape and the Type B target ai value are not from any
    % reference and are documented as such).
    %
    % IMPORTANT CAVEAT ON THESE EXPECTED VALUES: the typetxt/ai values
    % asserted below were NOT obtained by running ai_v2.m in MATLAB (this
    % assistant cannot execute MATLAB). They come from a from-scratch
    % Python reimplementation of ai_v2.m's exact algorithm (the fsg71
    % 7-point Savitzky-Golay derivative filter and the quantised
    % 4th-derivative zero-crossing search), built directly from this
    % file's source. That is weaker evidence than the hand-verified
    % values in testReadBPplusReaders.m, because the "oracle" here is a
    % reimplementation of the very function under test, not an
    % independent ground truth. Run this file in MATLAB and treat any
    % failure as informative either way - it may reveal a bug in the
    % Python reimplementation's understanding of fsg71/mminterp, not
    % necessarily a problem with ai_v2.m itself. Report back what
    % actually comes out.

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

        function testSyntheticTypeA(testCase)
            [~, metadata, ~, ~, ao] = read_BPplus(testCase.FixtureDir, ...
                'synthetic_ai_v2_TypeA.xml', testCase.Npoly, testCase.Frame);

            [ai, ~, ~, ~, ~, typetxt] = ai_v2(ao.p_av, metadata.samplerate);

            % Target ai/PP ratio (+0.29) matches Murgo 1980 Fig. 3's own
            % worked Type A example exactly.
            testCase.verifyEqual(typetxt, 'Type A');
            testCase.verifyEqual(ai, 29);
        end

        function testSyntheticTypeB(testCase)
            [~, metadata, ~, ~, ao] = read_BPplus(testCase.FixtureDir, ...
                'synthetic_ai_v2_TypeB.xml', testCase.Npoly, testCase.Frame);

            [ai, ~, ~, ~, ~, typetxt] = ai_v2(ao.p_av, metadata.samplerate);

            % No reference gives a numeric Type B example (Murgo 1980/1981
            % both describe it only as "intermediate...not illustrated",
            % bounded 0.0 < deltaP/PP < 0.12). ai=9 was chosen to sit
            % mid-band; it is not sourced from any reference.
            testCase.verifyEqual(typetxt, 'Type B');
            testCase.verifyEqual(ai, 9);
        end

        function testSyntheticTypeC(testCase)
            [~, metadata, ~, ~, ao] = read_BPplus(testCase.FixtureDir, ...
                'synthetic_ai_v2_TypeC.xml', testCase.Npoly, testCase.Frame);

            [ai, ~, ~, ~, ~, typetxt] = ai_v2(ao.p_av, metadata.samplerate);

            % Target ai/PP ratio (-0.19) matches Murgo 1980 Fig. 3's own
            % worked Type C example exactly.
            testCase.verifyEqual(typetxt, 'Type C');
            testCase.verifyEqual(ai, -19);
        end

    end
end
