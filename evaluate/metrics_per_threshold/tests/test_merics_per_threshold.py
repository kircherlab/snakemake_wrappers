from click.testing import CliRunner
import filecmp

from metrics_per_threshold import cli


def test1_metrics_per_threshold():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv', '--label-column', 'label',
                                 '--positive-label', 1, '--prediction-column', "prediction", '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest1.tsv")


def test2_metrics_per_threshold():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv', '--label-column', 'label',
                                 '--positive-label', 0, '--prediction-column', "prediction", '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest2.tsv")
