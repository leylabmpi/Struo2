#!/usr/bin/env python
from __future__ import print_function
import sys,os
import re
import argparse
import logging
import subprocess
import collections
from subprocess import Popen, PIPE


desc = 'Summarizing snakemake log'
epi = """DESCRIPTION:
Creates a table summarizing the snakemake run log.
"""
parser = argparse.ArgumentParser(description=desc,
                                 epilog=epi,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('log_file', metavar='log_file', type=str,
                    help='Snakemake screen log')
parser.add_argument('--acct', type=str, default='/var/lib/gridengine/default/common/accounting',
                    help='qacct accouting file (default: %(default)s)')
parser.add_argument('--version', action='version', version='0.0.1')

logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)

# regular expressions
## formatting
regex_esc = re.compile(r'\x1B\[[0-?]*[ -/]*[@-~]')
## search
regex_date = re.compile(r'^\[([A-Za-z]+ [A-Za-z]+ +\d+ \d\d:\d\d:\d\d \d{4})\]$')    
regex_rule = re.compile(r'^(localrule|rule) ([A-Za-z0-9-_]+):$')
regex_error_in = re.compile(r'^Error in rule .+:$')
regex_error_ex = re.compile(r'^Error executing rule .+see the cluster log$')
regex_finished = re.compile(r'^Finished job ([0-9]+)\.*$')
regex_exit = re.compile(r'Exiting because a job execution failed')
### search in rule section
regex_input = re.compile(r'^ +input: (.+)')
regex_log = re.compile(r'^ +log: (.+)')
regex_wildcards = re.compile(r'^ +wildcards: (.+)')
regex_resources = re.compile(r'^ +resources: (.+)')
regex_sm_job = re.compile(r'^ +jobid: ([0-9]+)')
regex_conda_env = re.compile(r'^ +conda-env: (.+)')
regex_clust_job = re.compile(r'^ *Submitted job [0-9]+ with external jobid \'([0-9]+)\'')
regex_error_clust_job = re.compile(r'^ +cluster_jobid: ([0-9]+)$')

def makehash():
    return collections.defaultdict(makehash)

def parse_rule(inF, D, rule_name, time_stamp):
    """Parsing single rule section
    """
    input_files = 'NA'
    log_file = 'NA'
    sm_job = 'NA'
    clust_job = 'NA'
    wildcards = 'NA'
    resources = 'NA'
    for line in inF:
        line = regex_esc.sub('', line.rstrip())
        if regex_date.search(line):
            # next section; adding to dict
            try:
                attempt_cnt = max(D[rule_name][input_files].keys()) + 1
            except (KeyError, ValueError) as e:
                attempt_cnt = 1
            X = {'log_file' : log_file,
                 'wildcards' : wildcards,
                 'resources' : resources,
                 'sm_job' : sm_job,
                 'clust_job' : clust_job,
                 'date_start' : time_stamp,
                 'date_end' : 'NA',
                 'pass_fail' : 'NA'}
            D[rule_name][input_files][attempt_cnt] = X
            return regex_date.search(line).group(1)
        elif regex_input.search(line):
            input_files = regex_input.search(line).group(1)
        elif regex_log.search(line):
            log_file = regex_log.search(line).group(1)
        elif regex_wildcards.search(line):
            wildcards = regex_wildcards.search(line).group(1)
            wildcards = wildcards.replace(' ', '')
        elif regex_resources.search(line):
            resources = regex_resources.search(line).group(1)
            resources = resources.replace(' ', '')
        elif regex_sm_job.search(line):
            sm_job = regex_sm_job.search(line).group(1)
        elif regex_clust_job.search(line):
            clust_job = regex_clust_job.search(line).group(1)

def parse_error(inF, D, rule_name, time_stamp):
    """Parsing error section
    """
    sm_job = 'NA'
    conda_env = 'NA'
    for line in inF:
        line = regex_esc.sub('', line.rstrip())
        if regex_date.search(line) or regex_exit.search(line):
            rule_name,input_files,attempt_cnt = find_sm_job(D, sm_job)
            D[rule_name][input_files][attempt_cnt]['pass_fail'] = 'error'
            D[rule_name][input_files][attempt_cnt]['date_end'] = time_stamp
            D[rule_name][input_files][attempt_cnt]['conda_env'] = conda_env
            try:
                return regex_date.search(line).group(1)
            except AttributeError:
                return None
        elif regex_sm_job.search(line):
            sm_job = regex_sm_job.search(line).group(1)
        elif regex_conda_env.search(line):
            conda_env = regex_conda_env.search(line).group(1)
            
def find_sm_job(D, jobID):
    """Getting snakemake job ID (if exists) from log file data structure
    """
    for rule_name,D1 in D.items():
        for input_files,D2 in D1.items():
            for attempt_cnt, D3 in sorted(D2.items(), reverse=True):
                try:
                    if D3['sm_job'] == jobID:
                        return rule_name, input_files, attempt_cnt
                except KeyError:
                    pass
    msg = 'WARNING: Cannot find rule for snakemake jobID: {}'
    logging.warning(msg.format(jobID))
    return 'NA', 'NA', 'NA'
            
def parse_section(inF, D, time_stamp):
    """Parsing section of log file (eg., rule or error)
    """
    if time_stamp is None:
        return None
    
    # next line
    line = inF.readline()
    line = regex_esc.sub('', line.rstrip())
    # what is next part?
    if line == '':
        pass
    elif regex_rule.search(line):
        rule_name = regex_rule.search(line).group(2)
        time_stamp = parse_rule(inF, D,
                                rule_name=rule_name,
                                time_stamp=time_stamp)
        parse_section(inF, D, time_stamp)
    elif regex_error_in.search(line):
        time_stamp = parse_error(inF, D, rule_name=line,
                                 time_stamp=time_stamp)
        parse_section(inF, D, time_stamp)
    elif regex_finished.search(line):
        sm_job = regex_finished.search(line).group(1)
        rule_name,input_files,attempt_cnt = find_sm_job(D, sm_job)
        D[rule_name][input_files][attempt_cnt]['pass_fail'] = 'finished'
        D[rule_name][input_files][attempt_cnt]['date_end'] = time_stamp
    elif regex_date.search(line):
        time_stamp = regex_date.search(line).group(1)
    else:
        msg = 'WARNING: Cannot parse: "{}"'
        logging.warning(msg.format(line))

def parse_log(log_file, D):
    """Parsing log file
    """
    logging.info('Parsing log: {}'.format(log_file))

    # parsing file
    with open(log_file) as inF:
        time_stamp = None
        rule_name = None
        attempt_cnt = None
        for line in inF:
            line = regex_esc.sub('', line.rstrip())
            # rule/finished/error
            if regex_date.search(line):
                time_stamp = regex_date.search(line).group(1)
                parse_section(inF, D, time_stamp)
    return D
            
def format_input_files(F):
    """Formatting 'input_files' field
    """
    F = F.split(', ')
    F = [os.path.split(x)[1] for x in F]
    return ','.join(F)[:40]

def load_qacct(acct_file, max_lines = 100000):
    """Loading tail of qacct into memory
    """
    logging.info('Loading qacct log: {}'.format(acct_file))
    p = Popen(['tac', acct_file], stdout=PIPE)
    output, err = p.communicate()
    rc = p.returncode

    qacct_info = {}
    for i,line in enumerate(output.decode().split('\n')):
        if i > max_lines:
            break
        line = line.split(':')
        try:
            # jobID : [exit_status, ru_wallclock, ru_utime, cpu, mem, io]
            qacct_info[line[5]] = [line[12], line[13], line[14], line[36], line[37], line[38]]
        except IndexError:
            qacct_info[line[5]] = ['', '', '', '', '', '']
    return qacct_info

def get_qacct_info(D, acct_file):
    """Getting job info from qacct log
    """
    # load acct info
    if os.path.isfile(acct_file):
        qacct_info = load_qacct(acct_file)
    else:
        logging.info('Cannot find qacct log: {}; Skipping'.format(acct_file))        
        return D

    # get exit code
    for rule_name,D1 in D.items():
        for input_files,D2 in D1.items():
            for attempt_cnt, D3 in D2.items():                
                try:
                    clust_job = D3['clust_job']
                except KeyError:
                    clust_job = None

                # info
                if clust_job is not None:
                    try:
                        x = qacct_info[clust_job]
                    except KeyError:
                        continue
                    except TypeError:
                        x = ['NA'] * 6
                    D[rule_name][input_files][attempt_cnt]['exit_status'] = x[0]
                    D[rule_name][input_files][attempt_cnt]['ru_wallclock'] = x[1]
                    D[rule_name][input_files][attempt_cnt]['ru_utime'] = x[2]
                    D[rule_name][input_files][attempt_cnt]['cpu'] = x[3]
                    D[rule_name][input_files][attempt_cnt]['mem'] = x[4]
                    D[rule_name][input_files][attempt_cnt]['io'] = x[5]

    return D        

def write_summary(D):
    """Writing summary table 
    """
    # header
    cols1 = ['finished_error', 'rule', 'input_files_trunc', 'attempt']
    cols2 = ['wildcards', 'sm_job', 'clust_job', 'exit_status', 'ru_wallclock', 'ru_utime',
             'cpu', 'mem', 'io', 'date_start', 'date_end', 'resources', 'log_file', 'conda_env']
    cols3 = ['input_files']
    print(','.join(cols1 + cols2 + cols3))
    
    # body
    for rule_name,D1 in D.items():
        for input_files,D2 in D1.items():
            for attempt_cnt, D3 in D2.items():
                # pass fail
                try:
                    line = [D3['pass_fail']]
                except KeyError:
                    line = ['NA']
                # main keys
                line += [rule_name,
                         format_input_files(input_files),
                         attempt_cnt]
                # extra info
                for x in cols2:
                    try:
                        line.append(D3[x])
                    except KeyError:
                        line.append('NA')
                line.append(input_files[:5000].replace(' ', ''))
                # writing
                print(','.join([str(x).replace(',', ':') for x in line]))

def main(args):
    """Main interface
    """
    # create log file data structure
    D = makehash()
    # parsing snakemake log file
    D = parse_log(args.log_file, D)
    # getting info from qacct log table
    D = get_qacct_info(D, args.acct)
    # writing summary table
    write_summary(D)


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)
