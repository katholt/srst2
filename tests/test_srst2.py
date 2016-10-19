#!/usr/bin/env python

import os
import sys
import unittest

from mock import MagicMock, patch
from StringIO import StringIO

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'scripts')))

import srst2

class TestGetSamtoolsExec(unittest.TestCase):
	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_samtools_with_overide(self, env_mock, path_mock):
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_SAMTOOLS': '/usr/bin/samtools'}
		env_mock.get.side_effect = fake_env_variables.get
		samtools_exec = srst2.get_samtools_exec()
		self.assertEqual(samtools_exec, '/usr/bin/samtools')

	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_samtools_with_default(self, env_mock, path_mock):
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {}
		env_mock.get.side_effect = fake_env_variables.get
		samtools_exec = srst2.get_samtools_exec()
		self.assertEqual(samtools_exec, 'samtools')

	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_samtools_with_missing(self, env_mock, path_mock):
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_SAMTOOLS': '/missing/samtools'}
		env_mock.get.side_effect = fake_env_variables.get
		samtools_exec = srst2.get_samtools_exec()
		self.assertEqual(samtools_exec, 'samtools')

class TestGetBowtieExecs(unittest.TestCase):
	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_bowtie_with_overides(self, env_mock, path_mock):
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2': '/usr/bin/bowtie2',
							  'SRST2_BOWTIE2_BUILD': '/usr/bin/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		bowtie_exec, bowtie_build_exec = srst2.get_bowtie_execs()
		self.assertEqual(bowtie_exec, '/usr/bin/bowtie2')
		self.assertEqual(bowtie_build_exec, '/usr/bin/bowtie2-build')

	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_bowtie_with_defaults(self, env_mock, path_mock):
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {}
		env_mock.get.side_effect = fake_env_variables.get
		bowtie_exec, bowtie_build_exec = srst2.get_bowtie_execs()
		self.assertEqual(bowtie_exec, 'bowtie2')
		self.assertEqual(bowtie_build_exec, 'bowtie2-build')

	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_bowtie_with_mixture(self, env_mock, path_mock):
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2': '/usr/bin/bowtie2'}
		env_mock.get.side_effect = fake_env_variables.get
		bowtie_exec, bowtie_build_exec = srst2.get_bowtie_execs()
		self.assertEqual(bowtie_exec, '/usr/bin/bowtie2')
		self.assertEqual(bowtie_build_exec, '/usr/bin/bowtie2-build')

	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_bowtie_with_other_mixture(self, env_mock, path_mock):
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2_BUILD': '/usr/bin/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		bowtie_exec, bowtie_build_exec = srst2.get_bowtie_execs()
		self.assertEqual(bowtie_exec, 'bowtie2')
		self.assertEqual(bowtie_build_exec, '/usr/bin/bowtie2-build')

	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_bowtie_with_missing(self, env_mock, path_mock):
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2': '/missing/bowtie2',
							  'SRST2_BOWTIE2_BUILD': '/usr/bin/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		bowtie_exec, bowtie_build_exec = srst2.get_bowtie_execs()
		self.assertEqual(bowtie_exec, 'bowtie2')
		self.assertEqual(bowtie_build_exec, '/usr/bin/bowtie2-build')

	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_bowtie_with_other_missing(self, env_mock, path_mock):
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2': '/usr/bin/bowtie2',
							  'SRST2_BOWTIE2_BUILD': '/missing/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		bowtie_exec, bowtie_build_exec = srst2.get_bowtie_execs()
		self.assertEqual(bowtie_exec, '/usr/bin/bowtie2')
		self.assertEqual(bowtie_build_exec, '/usr/bin/bowtie2-build')

class TestBowtieIndex(unittest.TestCase):
	@patch('srst2.logging')
	@patch('srst2.check_command_versions')
	@patch('srst2.run_command')
	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_bowtie_index_with_overides(self, env_mock, path_mock, run_mock,
										version_mock, logging_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {'SRST2_BOWTIE2': '/usr/bin/bowtie2',
							  'SRST2_BOWTIE2_BUILD': '/usr/bin/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		srst2.bowtie_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args_list[0][0][0], ['/usr/bin/bowtie2',
														'--version'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['/usr/bin/bowtie2-build', 'foo', 'foo'])

	@patch('srst2.logging')
	@patch('srst2.check_command_versions')
	@patch('srst2.run_command')
	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_bowtie_index_with_defaults(self, env_mock, path_mock, run_mock,
										version_mock, logging_mock):
		path_mock.exists.return_value = False
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fake_env_variables = {}
		env_mock.get.side_effect = fake_env_variables.get
		srst2.bowtie_index(['foo'])
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args_list[0][0][0], ['bowtie2',
														'--version'])
		self.assertEqual(run_mock.call_count, 1)
		run_mock.assert_called_once_with(['bowtie2-build', 'foo', 'foo'])

class TestRunBowtie(unittest.TestCase):
	@patch('srst2.logging')
	@patch('srst2.check_command_versions')
	@patch('srst2.run_command')
	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_run_bowtie_with_overide(self, env_mock, path_mock, run_mock,
									 version_mock, logging_mock):
		fake_env_variables = {'SRST2_SAMTOOLS': '/usr/bin/samtools',
							  'SRST2_BOWTIE2': '/usr/bin/bowtie2',
							  'SRST2_BOWTIE2_BUILD': '/usr/bin/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		arg_mock = MagicMock()
		arg_mock.read_type = 'foo'
		arg_mock.stop_after = False
		arg_mock.other = False
		arg_mock.threads = 4
		arg_mock.use_existing_bowtie2_sam = False
		actual_sam = srst2.run_bowtie('mapping_file', 'sample', ['fastq'], arg_mock,
									  'db_name', 'db_path')
		self.assertEqual(actual_sam, 'mapping_file.sam')
		self.assertEqual(version_mock.call_count, 2)
		self.assertEqual(version_mock.call_args_list[0][0][0], ['/usr/bin/bowtie2',
														'--version'])
		self.assertEqual(version_mock.call_args_list[1][0][0], ['/usr/bin/samtools'])
		expected_bowtie2_command = [
			'/usr/bin/bowtie2',
			'-U', 'fastq',
			'-S', 'mapping_file.sam',
			'-foo',
			'--very-sensitive-local',
			'--no-unal',
			'-a',
			'-x', 'db_path',
			'--threads', '4'
		]
		run_mock.assert_called_once_with(expected_bowtie2_command)

	@patch('srst2.logging')
	@patch('srst2.check_command_versions')
	@patch('srst2.run_command')
	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_run_bowtie_with_defaults(self, env_mock, path_mock, run_mock,
									  version_mock, logging_mock):
		fake_env_variables = {}
		env_mock.get.side_effect = fake_env_variables.get
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		arg_mock = MagicMock()
		arg_mock.read_type = 'foo'
		arg_mock.stop_after = False
		arg_mock.other = False
		arg_mock.threads = 4
		arg_mock.use_existing_bowtie2_sam = False
		actual_sam = srst2.run_bowtie('mapping_file', 'sample', ['fastq'], arg_mock,
									  'db_name', 'db_path')
		self.assertEqual(actual_sam, 'mapping_file.sam')
		self.assertEqual(version_mock.call_count, 2)
		self.assertEqual(version_mock.call_args_list[0][0][0], ['bowtie2',
														'--version'])
		self.assertEqual(version_mock.call_args_list[1][0][0], ['samtools'])
		expected_bowtie2_command = [
			'bowtie2',
			'-U', 'fastq',
			'-S', 'mapping_file.sam',
			'-foo',
			'--very-sensitive-local',
			'--no-unal',
			'-a',
			'-x', 'db_path',
			'--threads', '4'
		]
		run_mock.assert_called_once_with(expected_bowtie2_command)

class TestMPileup(unittest.TestCase):
	@patch('srst2.open', create=True)
	@patch('srst2.logging')
	@patch('srst2.check_command_versions')
	@patch('srst2.run_command')
	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_get_pileup_with_overides(self, env_mock, path_mock, run_mock,
									  version_mock, logging_mock, open_mock):
		fake_env_variables = {'SRST2_SAMTOOLS': '/usr/bin/samtools',
							  'SRST2_BOWTIE2': '/usr/bin/bowtie2',
							  'SRST2_BOWTIE2_BUILD': '/usr/bin/bowtie2-build'}
		env_mock.get.side_effect = fake_env_variables.get
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		arg_mock = MagicMock()
		arg_mock.mapq = 30
		arg_mock.baseq = 40
		arg_mock.samtools_args = []
		arg_mock.keep_interim_alignment = True # They're not actually created
		fake_file = MagicMock()
		fake_open_context = MagicMock(**{'__enter__.return_value': fake_file})
		open_mock.return_value = fake_open_context
		srst2.get_pileup(arg_mock, 'mapping_file', 'raw_bowtie_sam',
						 'bowtie_sam_mod', 'fasta', 'pileup')
		
		expected_samtools_view_command = [
			'/usr/bin/samtools',
			'view',
			'-b',
			'-o', 'mapping_file.unsorted.bam',
			'-q', '30',
			'-S', 'bowtie_sam_mod'
		]
		run_mock.assert_any_call(expected_samtools_view_command)
		
		expected_samtools_sort_command = [
			'/usr/bin/samtools',
			'sort',
			'mapping_file.unsorted.bam',
			'mapping_file.sorted'
		]
		run_mock.assert_any_call(expected_samtools_sort_command)

		expected_mpileup_command = [
			'/usr/bin/samtools',
			'mpileup',
			'-L', '1000',
			'-f', 'fasta',
			'-Q', '40',
			'-q', '30',
			'-B', 'mapping_file.sorted.bam'
		]
		run_mock.assert_any_call(expected_mpileup_command,
									stdout=fake_file)

	@patch('srst2.open', create=True)
	@patch('srst2.logging')
	@patch('srst2.check_command_versions')
	@patch('srst2.run_command')
	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_get_pileup_with_defaults(self, env_mock, path_mock, run_mock,
									  version_mock, logging_mock, open_mock):
		fake_env_variables = {}
		env_mock.get.side_effect = fake_env_variables.get
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		arg_mock = MagicMock()
		arg_mock.mapq = 30
		arg_mock.baseq = 40
		arg_mock.samtools_args = []
		arg_mock.keep_interim_alignment = True # They're not actually created
		fake_file = MagicMock()
		fake_open_context = MagicMock(**{'__enter__.return_value': fake_file})
		open_mock.return_value = fake_open_context
		srst2.get_pileup(arg_mock, 'mapping_file', 'raw_bowtie_sam',
						 'bowtie_sam_mod', 'fasta', 'pileup')
		
		expected_samtools_view_command = [
			'samtools',
			'view',
			'-b',
			'-o', 'mapping_file.unsorted.bam',
			'-q', '30',
			'-S', 'bowtie_sam_mod'
		]
		run_mock.assert_any_call(expected_samtools_view_command)
		
		expected_samtools_sort_command = [
			'samtools',
			'sort',
			'mapping_file.unsorted.bam',
			'mapping_file.sorted'
		]
		run_mock.assert_any_call(expected_samtools_sort_command)

		expected_mpileup_command = [
			'samtools',
			'mpileup',
			'-L', '1000',
			'-f', 'fasta',
			'-Q', '40',
			'-q', '30',
			'-B', 'mapping_file.sorted.bam'
		]
		run_mock.assert_any_call(expected_mpileup_command,
									stdout=fake_file)

class TestSamtoolsIndex(unittest.TestCase):
	@patch('srst2.check_command_versions')
	@patch('srst2.run_command')
	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_samtools_index_with_overides(self, env_mock, path_mock, run_mock,
									   version_mock):
		fake_env_variables = {'SRST2_SAMTOOLS': '/usr/bin/samtools'}
		path_mock.exists.return_value = False
		env_mock.get.side_effect = fake_env_variables.get
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fai_file = srst2.samtools_index('fasta')
		
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args_list[0][0][0], ['/usr/bin/samtools'])
		expected_samtools_command = [
			'/usr/bin/samtools',
			'faidx',
			'fasta'
		]
		run_mock.assert_called_once_with(expected_samtools_command)

	@patch('srst2.check_command_versions')
	@patch('srst2.run_command')
	@patch('srst2.os.path')
	@patch('srst2.os.environ')
	def test_samtools_index_with_defaults(self, env_mock, path_mock, run_mock,
									   version_mock):
		fake_env_variables = {}
		path_mock.exists.return_value = False
		env_mock.get.side_effect = fake_env_variables.get
		path_mock.isfile.side_effect = lambda f: 'missing' not in f
		fai_file = srst2.samtools_index('fasta')
		
		self.assertEqual(version_mock.call_count, 1)
		self.assertEqual(version_mock.call_args_list[0][0][0], ['samtools'])
		expected_samtools_command = [
			'samtools',
			'faidx',
			'fasta'
		]
		run_mock.assert_called_once_with(expected_samtools_command)

if __name__ == '__main__':
	unittest.main()
