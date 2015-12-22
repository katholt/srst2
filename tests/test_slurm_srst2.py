#!/usr/bin/env python

import os
import sys
import unittest

from mock import patch

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'scripts')))

import slurm_srst2

class TestSamtoolExec(unittest.TestCase):
	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_samtools_index_with_overide(self, env_mock, path_mock, run_mock,
										 version_mock, stdout_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		env_mock.get.side_effect = {'SRST2_SAMTOOLS': '/usr/bin/samtools'}.get
		slurm_srst2.samtools_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['/usr/bin/samtools'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['/usr/bin/samtools', 'faidx',
										  'foo'])

	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_samtools_index_with_default(self, env_mock, path_mock, run_mock,
										 version_mock, stdout_mock):
		path_mock.exists.return_value = False
		env_mock.get.return_value = None
		slurm_srst2.samtools_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['samtools'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['samtools', 'faidx',
										  'foo'])

	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_samtools_index_with_missing_overide(self, env_mock, path_mock,
												 run_mock, version_mock,
												 stdout_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		env_mock.get.side_effect = {'SRST2_SAMTOOLS': '/missing/samtools'}.get
		slurm_srst2.samtools_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['samtools'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['samtools', 'faidx',
										  'foo'])
class TestBowtieIndex(unittest.TestCase):
	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_bowtie_index_with_overides(self, env_mock, path_mock, run_mock,
										version_mock, stdout_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2': '/usr/bin/bowtie2',
							  'SRST2_BOWTIE2_BUILD': '/usr/bin/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		slurm_srst2.bowtie_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['/usr/bin/bowtie2',
														'--version'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['/usr/bin/bowtie2-build', 'foo', 'foo'])

	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_bowtie_index_with_partial_overide(self, env_mock, path_mock, run_mock,
											   version_mock, stdout_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2': '/usr/bin/bowtie2'}
		env_mock.get.side_effect = fake_env_variables.get
		slurm_srst2.bowtie_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['/usr/bin/bowtie2',
														'--version'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['/usr/bin/bowtie2-build', 'foo', 'foo'])

	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_bowtie_index_with_other_overide(self, env_mock, path_mock, run_mock,
											 version_mock, stdout_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2_BUILD': '/usr/bin/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		slurm_srst2.bowtie_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['bowtie2',
														'--version'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['/usr/bin/bowtie2-build', 'foo', 'foo'])

	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_bowtie_index_with_default(self, env_mock, path_mock, run_mock,
										 version_mock, stdout_mock):
		path_mock.exists.return_value = False
		env_mock.get.return_value = None
		slurm_srst2.bowtie_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['bowtie2',
														'--version'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['bowtie2-build', 'foo', 'foo'])

	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_bowtie_index_with_missing_overide(self, env_mock, path_mock,
												 run_mock, version_mock,
												 stdout_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		env_mock.get.side_effect = {'SRST2_SAMTOOLS': '/missing/bowtie'}.get
		slurm_srst2.bowtie_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['bowtie2',
														'--version'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['bowtie2-build', 'foo', 'foo'])

	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_bowtie_index_with_one_missing_overide(self, env_mock, path_mock,
												   run_mock, version_mock,
												   stdout_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2': '/missing/bowtie2',
							  'SRST2_BOWTIE2_BUILD': '/usr/bin/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		slurm_srst2.bowtie_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['bowtie2',
														'--version'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['/usr/bin/bowtie2-build', 'foo', 'foo'])

	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_bowtie_index_with_other_missing_overide(self, env_mock, path_mock,
													 run_mock, version_mock,
													 stdout_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2': '/usr/bin/bowtie2',
							  'SRST2_BOWTIE2_BUILD': '/missing/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		slurm_srst2.bowtie_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['/usr/bin/bowtie2',
														'--version'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['/usr/bin/bowtie2-build', 'foo', 'foo'])

	@patch('slurm_srst2.sys.stdout')
	@patch('slurm_srst2.check_command_versions')
	@patch('slurm_srst2.run_command')
	@patch('slurm_srst2.os.path')
	@patch('slurm_srst2.os.environ')
	def test_bowtie_index_with_missing_inference(self, env_mock, path_mock,
												 run_mock, version_mock,
												 stdout_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: f == '/usr/bin/bowtie2'
		fake_env_variables = {'SRST2_BOWTIE2': '/usr/bin/bowtie2'}
		env_mock.get.side_effect = fake_env_variables.get
		slurm_srst2.bowtie_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args[0][0], ['/usr/bin/bowtie2',
														'--version'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['bowtie2-build', 'foo', 'foo'])

if __name__ == '__main__':
	unittest.main()
