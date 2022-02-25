from click.testing import CliRunner
import filecmp

from metrics import cli


def test1_auc():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv', '--label-column', 'label',
                                 '--positive-label', 1, '--prediction-column', "prediction", '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest1.tsv")


def test2_auc():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv', '--label-column', 'label',
                                 '--positive-label', 0, '--prediction-column', "prediction", '--output', "tests/output.tsv"])
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest2.tsv")
