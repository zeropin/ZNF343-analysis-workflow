#!/usr/local/bin/python3

import sys
import re
import math

def readBED(filename):
        BED_file = open(filename, 'r');
        BED = {};
        key = '';

        for line in BED_file:
                element = line.split();

                if element[5] == '+':
                        key = element[0]+'+'+element[1];
                elif element[5] == '-':
                        key = element[0]+'-'+element[2];

                if key in BED:
                        BED[key] += 1;
                else:
                        BED[key] = 1;

        BED_file.close();
        return BED;

def detect(anchors, reads, anchor, border_width=40, forward_shift=1, backward_shift=30, output_file='ZNF343.exo'): 
	aligned_peaks = open(anchors,'r');
	bed	= readBED(reads);
	output  = open(output_file, 'w');

	anchor_length = len(anchor);

	seq_count = 0;
	direction = '';

	len_total = 2*border_width+anchor_length;

	Forward_count = [0]*(len_total);
	Reverse_count = [0]*(len_total);

	for line in aligned_peaks:
		if line[0]=='>':
			words = re.split('>|:|-|\(|\)',line);
			chromosome = words[1];
			start = int(words[2])+int(words[4]);
			end = int(words[2])+int(words[5]);
		
			if words[6]=='+':
				direction = '+';
			elif words[6]=='':
				direction = '-';
		else:
			sequence = line;
			seq_count += 1;

			forward_reads = 0;
			reverse_reads = 0;

			if direction=='+':
				for i in range(start-border_width, end+border_width+1):
					key = chromosome+'+'+str(i);
					key2 = chromosome+'-'+str(i);
					if (key in bed):
						Forward_count[i-start+border_width] += bed[key];
						if i<end:
							forward_reads += bed[key];
					if (key2 in bed):
						Reverse_count[i-start+border_width] += bed[key2];
						if i>start:
							reverse_reads += bed[key2];
				exo_reads = forward_reads + reverse_reads;
				output.write(chromosome+'	'+str(start-backward_shift-1)+'	'+str(end+forward_shift)+'	'+str(forward_reads+reverse_reads)+'	'+'name'+'	'+direction+'\n');

			elif direction=='-':
				for i in range(end+border_width,start-border_width-1,-1):
					key = chromosome+'-'+str(i);
					key2 = chromosome+'+'+str(i);
					if (key in bed):
						Forward_count[end+border_width-i] += bed[key];
						if i>start:
							forward_reads += bed[key];
					if (key2 in bed):
						Reverse_count[end+border_width-i] += bed[key2];
						if i<end:
							reverse_reads += bed[key2];
				output.write(chromosome+'	'+str(start-forward_shift-1)+'	'+str(end+backward_shift)+'	'+str(forward_reads+reverse_reads)+'	'+'name'+'	'+direction+'\n');
			
							
	output.close();

	background_forward = sum(Forward_count[len_total-6:len_total])*1.0/6;
	background_reverse = sum(Reverse_count[0:6])*1.0/6;

	subtraction_forward = background_forward*(border_width+anchor_length-1)/seq_count;
	subtraction_reverse = background_reverse*(border_width+anchor_length-1)/seq_count;

	print("Total number of sequences:"+str(seq_count));
	print("Length of tracked region:"+str(len_total));
	print("Average background forward reads per position:"+'%.1f'%background_forward);
	print("Average background reverse reads per position:"+'%.1f'%background_reverse);
	print("Background value to be subtracted for forward/reverse reads per sequence:"+'%.1f'%subtraction_forward+'	'+'%.1f'%subtraction_reverse);
	
	return [Forward_count, Reverse_count];

