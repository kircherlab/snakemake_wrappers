from click.testing import CliRunner
import filecmp

from summarize_columns import cli


def test1_summarize_columns():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A', '--column', 'B',
                                 '--new-column-name', 'average',
                                 '--operation', "mean", '--operation', "abs_mean", 
                                 '--output', "tests/output.tsv"])                          
    assert result.exit_code == 1

def test2_summarize_columns():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A', '--column', 'B',
                                 '--new-column-name', 'average',
                                 '--operation', "abs", 
                                 '--output', "tests/output.tsv"])                          
    assert result.exit_code == 1


def test3_summarize_columns():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A',
                                 '--new-column-name', 'max',
                                 '--operation', "max",
                                 '--output', "tests/output.tsv"])                      
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest3.tsv")

def test4_summarize_columns():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A', '--column', 'B',
                                 '--new-column-name', 'max',
                                 '--operation', "max",
                                 '--output', "tests/output.tsv"])                      
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest4.tsv")

def test5_summarize_columns():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A', '--column', 'B',
                                 '--new-column-name', 'mean',
                                 '--operation', "mean",
                                 '--output', "tests/output.tsv"])                      
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest5.tsv")

def test6_summarize_columns():
    runner = CliRunner()
    result = runner.invoke(cli, ['--input', 'tests/input.tsv',
                                 '--column', 'A', '--column', 'B',
                                 '--new-column-name', 'abs_mean',
                                 '--operation', "abs_mean",
                                 '--output', "tests/output.tsv"])                      
    assert result.exit_code == 0
    assert filecmp.cmp("tests/output.tsv", "tests/outTest6.tsv")