# Last updated Jan 2025
# Authors: Auden Cote-L'Heureux, Mario Ceron-Romero.

# This script contains the entirety of the contamination loop, an iterative tool to assess
# and remove contamination based on analyses of single gene trees. This tool allows for the
# use of one or both of two main methods of contamination assessment informed by tree topology.
# The first method – “clade-based contamination removal” –  is intended for cases when the
# user is interested in genes present in a group of organisms with multiple representative
# samples and/or species in the gene trees. For a given set of target taxa, this method identifies
# robust, monophyletic clades containing those taxa within each gene tree, and removes all other
# sequences belonging to the group of target taxa. The second available method is 'sister-based
# contamination removal', where sequences are removed based on their sister relationship to
# user-specified contaminants. After sequences are removed in each iteration, the loop re-aligns and
# re-builds each tree excluding all removed sequences, and then repeats a specified number of times
# (--n_loops). Running in botb of these modes requires the user to input several parameters,
# which can be found in the manual, and in the 'CL' argument group in utils.py.

#Dependencies
import os, sys, re
from Bio import SeqIO
import ete3
import guidance
import trees
from statistics import mean

def get_newick(fname):
	"""Extract the Newick string from a Nexus file.

	Args:
		fname (str): The path to the Nexus file.

	Returns:
		str: The Newick string extracted from the file, or an empty string if no Newick string is found.
	"""

	newick = ''
	for line in open(fname):
		line = line.split(' ')[-1]
		if(line.startswith('(') or line.startswith('tree1=')):
			newick = line.split('tree1=')[-1].replace("'", '').replace('\\', '')

	return newick

def reroot(tree):
	"""Reroot the tree on the largest clade of a specified taxonomic group.

	Search taxonomic groups in the following order: Ba/Za, Op, Pl, Am, Ex, Sr.
	Utilize the first group with a qualifying clade as the outgroup. Leave the tree unchanged if none qualifies.

	Args:
		tree (ete3.TreeNode): The input tree to be rerooted.

	Returns:
		ete3.TreeNode: The rerooted tree, or the original tree if no suitable clade is found.
	"""
	def get_best_clade(taxon):
		"""Find the largest clade of a given taxonomic group.

		Examine nodes with at least three leaves, excluding those that are subsumed by previously accepted clades.
		A qualifying clade is one that contains the maximum number of leaves from the specified taxonomic group,
		while having at most two leaves from other groups.

		Args:
			taxon (tuple): A tuple of taxonomic group codes (e.g., ('Ba', 'Za')).

		Returns:
			ete3.TreeNode: The largest qualifying clade of the specified taxonomic group, or an empty list if none is found.
			Note: The clade is represented as a TreeNode object, not a list of leaf names.
		"""
		best_size = 0; best_clade = []; seen_leaves = []
		#Traverse all nodes
		for node in tree.traverse('levelorder'):
			#If the node is big enough and not subsumed by a node we've already accepted
			if len(node) >= 3 and len(list(set(seen_leaves) & set([leaf.name for leaf in node]))) == 0:
				leaves = [leaf.name for leaf in node]

				#Create a record of leaves that belong to the taxonomic group
				target_leaves = set()
				for leaf in leaves[::-1]:
					if leaf[:2] in taxon:
						target_leaves.add(leaf[:10])
						leaves.remove(leaf)

				#If this clade is better than any clade we've seen before, grab it
				if len(target_leaves) > best_size and len(leaves) <= 2:
					best_clade = node
					best_size = len(target_leaves)
					seen_leaves.extend([leaf.name for leaf in node])

		return best_clade

	#Get the biggest clade for each taxonomic group (stops once it finds one)
	for taxon in [('Ba', 'Za'), ('Op'), ('Pl'), ('Am'), ('Ex'), ('Sr')]:
		clade = get_best_clade(taxon)
		if len([leaf for leaf in clade if leaf.name[:2] in taxon]) > 3:
			tree.set_outgroup( clade)

			break

	return tree

def get_subtrees(args, file):
	"""Get the sequences to keep from a tree based on clade grabbing rules.

	Read and reroot a gene tree, then select clades that meet the specified taxonomic requirements
	and contamination thresholds. A clade is selected if 1) each of its children contains at least
	one target taxon; 2) it has at least min_target_presence distinct target taxa; 3) it has at least
	required_taxa_num distinct required taxa; and 4) its number of non-target leaves is at most
	num_contams (or, if num_contams < 1, at most that fraction of the number of target taxa).
	Avoid selecting overlapping clades, and retain non-target sequences and those listed as
	exceptions. Rules come from args.clade_grabbing_rules_file if given, and otherwise from the
	individual command-line arguments. Exit the program if a rules or exceptions file cannot be found.

	Args:
		args (argparse.Namespace): Command-line arguments containing clade grabbing rules and exceptions.
		file (str): Path to the input tree file.

	Returns:
		list[str]: Names of sequences to keep based on the clade grabbing rules.
	"""
	newick = get_newick(file)

	tree = ete3.Tree(newick)

	try:
		tree = reroot(tree)
	except:
		print('\nUnable to re-root the tree ' + file + ' (maybe it had only 1 major clade, or an inconvenient polytomy). Skipping this step and continuing to try to grab robust clades from the tree.\n')

	#Reading in clade grabbing exceptions
	exceptions = []
	if args.clade_grabbing_exceptions != None:
		if os.path.isfile(args.clade_grabbing_exceptions):
			exceptions = [line.strip() for line in open(args.clade_grabbing_exceptions)]
		else:
			print('\nError: it looks like you tried to input a clade grabbing exceptions file, but it could not be found.\n')
			exit()

	#Reading clade grabbing rules
	rules_per_clade = []
	if args.clade_grabbing_rules_file != None:
		if os.path.isfile(args.clade_grabbing_rules_file):
			lines = [line.strip().split('\t') for line in open(args.clade_grabbing_rules_file) if len(line.strip().split('\t')) == 5]

			for line in lines:
				if line[4].lower() == 'na':
					rules_per_clade.append({ 'target_taxa' : line[0], 'num_contams' : float(line[1]), 'min_target_presence' : int(line[2]), 'required_taxa' : line[3], 'required_taxa_num' : 0 })
				else:
					rules_per_clade.append({ 'target_taxa' : line[0], 'num_contams' : float(line[1]), 'min_target_presence' : int(line[2]), 'required_taxa' : line[3], 'required_taxa_num' : int(line[4]) })

		else:
			print('\nError: it looks like you tried to input a clade grabbing rules file, but it could not be found.\n')
			exit()
	else:
		rules_per_clade.append({ 'target_taxa' : args.target_taxa, 'num_contams' : args.num_contams, 'min_target_presence' : args.min_target_presence, 'required_taxa' : args.required_taxa, 'required_taxa_num' : args.required_taxa_num })

	#Reformatting rules
	for clade in rules_per_clade:
		if os.path.isfile(clade['target_taxa']):
			try:
				clade['target_taxa'] = [l.strip() for l in open(clade['target_taxa']).readlines() if l.strip() != '']
			except AttributeError:
				print('\nError: invalid "target_taxa" input (' + clade['target_taxa'] + '). This must be a comma-separated list of any number of digits/characters to describe focal taxa (e.g. Sr_ci_S,Am_t), or a file with the extension .txt containing a list of complete or partial taxon codes. All sequences containing the complete/partial code will be identified as belonging to target taxa.\n')
				exit()
		else:
			clade['target_taxa'] = [code.strip() for code in clade['target_taxa'].split(',') if code.strip() != '']

		if clade['required_taxa'] != None:
			if os.path.isfile(clade['required_taxa']):
				try:
					clade['required_taxa'] = [l.strip() for l in open(clade['required_taxa']).readlines() if l.strip() != '']
				except AttributeError:
					print('\nError: invalid "required_taxa" argument. This must be a comma-separated list of any number of digits/characters (e.g. Sr_ci_S,Am_t), or a file with the extension .txt containing a list of complete or partial taxon codes, to describe taxa that MUST be present in a clade for it to be selected (e.g. you may want at least one whole genome).\n')
			else:
				clade['required_taxa'] = [code.strip() for code in clade['required_taxa'].split(',') if code.strip() != '' and code.strip().lower() != 'na']

			clade['target_taxa'] = clade['target_taxa'] + clade['required_taxa']
		else:
			clade['required_taxa'] = []

	#Creating a record of selected subtrees, and all of the leaves in those subtrees
	selected_leaves = []

	#For each set of rules (set of target taxa)
	for clade in rules_per_clade:
		seen_leaves = []

		#Iterating through all nodes in tree, starting at "root" then working towards leaves
		for node in tree.traverse('levelorder'):
			#If a node is large enough and is not contained in an already selected clade
			if len(node) >= clade['min_target_presence'] and len(list(set(seen_leaves) & set([leaf.name for leaf in node]))) == 0:
				leaves = [leaf.name for leaf in node]

				#Accounting for cases where e.g. one child is a contaminant, and the other child is a good clade with 1 fewer than the max number of contaminants
				children_keep = 0
				for child in node.children:
					for code in clade['target_taxa']:
						taken = False
						for leaf in child:
							if leaf.name.startswith(code):
								children_keep += 1
								taken = True
								break

						if taken:
							break

				if children_keep == len(node.children):

					#Creating a record of all leaves belonging to the target/"at least" group of taxa, and any other leaves are contaminants
					target_leaves = set(); at_least_leaves = set(); target_leaves_full_names = []
					for code in clade['target_taxa']:
						for leaf in leaves[::-1]:
							if leaf.startswith(code):
								target_leaves.add(leaf[:10])
								target_leaves_full_names.append(leaf)

								for req in clade['required_taxa']:
									if leaf.startswith(req):
										at_least_leaves.add(leaf[:10])
										break

								leaves.remove(leaf)

					#Grab a clade as a subtree if 1) it has enough target taxa; 2) it has enough "at least" taxa; 3) it does not have too many contaminants
					if len(target_leaves) >= clade['min_target_presence'] and len(at_least_leaves) >= clade['required_taxa_num'] and ((clade['num_contams'] < 1 and len(leaves) <= clade['num_contams'] * len(target_leaves)) or len(leaves) <= clade['num_contams']):
						selected_leaves.extend(target_leaves_full_names)

						seen_leaves.extend([leaf.name for leaf in node])

	all_clades = [clade for group in rules_per_clade for clade in group['target_taxa']]

	seqs2keep = [leaf.name for leaf in tree if leaf.name in selected_leaves or not any([leaf.name.startswith(clade) for clade in all_clades]) or any([leaf.name.startswith(ex) for ex in exceptions])]

	return seqs2keep

def get_sisters(args, file, sister_contam_per_tax, subsister_contam_per_tax):
	"""Identify sequences to keep based on sister and sub-sister relationships.

	Read and reroot a gene tree. For each sequence matching a taxon in the rules, move towards
	the root until reaching a node that contains a leaf from a different taxon (taxa grouped
	together in args.cocontaminants are treated as one); the other taxa in that node are the
	sisters. The sub-sisters are the leaves in the next node up that are not in the sister node
	(only if args.subsister_rules is given). Remove a sequence if 1) all of its sisters are listed
	as contaminants for its taxon and its branch length is at most the smallest applicable
	threshold times the mean terminal branch length of the tree; or 2) all of its sub-sisters are listed
	as contaminants for its taxon. Exit the program if the co-contaminants file cannot be found.

	Args:
		args (argparse.Namespace): Command-line arguments for contamination removal, including the co-contaminants file and sub-sister rules file.
		file (str): Path to the input tree file.
		sister_contam_per_tax (dict): Sister contamination rules for each taxon.
		subsister_contam_per_tax (dict): Sub-sister contamination rules for each taxon.

	Returns:
		list[str]: Names of sequences to keep.
	"""

	seqs2remove = []

	#Read the tree using ete3 and reroot it using the above function
	newick = get_newick(file)
	tree = ete3.Tree(newick)

	try:
		tree = reroot(tree)
	except:
		print('\nUnable to re-root the tree ' + file + ' (maybe it had only 1 major clade, or an inconvenient polytomy). Skipping this step and continuing to try to grab robust clades from the tree.\n')

	mean_bl = mean([leaf.dist for leaf in tree])

	all_taxa_in_tree = [leaf.name[:10] for leaf in tree]

	#Reading in cocontaminants
	coconts = { }
	if args.cocontaminants != None:
		if os.path.isfile(args.cocontaminants):
			for line in open(args.cocontaminants):
				line = line.strip().split('\t')
				if len(line) == 2:
					if line[0] in all_taxa_in_tree:
						coconts.update({ line[0] : line[1] })
		else:
			print('\nERROR: It looks like you tried to input a co-contaminants file to the contamination loop, but the file could not be found.\n')
			exit()

	for tax in all_taxa_in_tree:
		if tax not in coconts:
			coconts.update({ tax : tax })

	#For each sequence
	for leaf in tree:
		bad_sisters = { contam[0] : contam[1] for tax in sister_contam_per_tax for contam in sister_contam_per_tax[tax] if leaf.name.startswith(tax) }
		bad_subsisters = [contam for tax in subsister_contam_per_tax for contam in subsister_contam_per_tax[tax] if leaf.name.startswith(tax)]

		if len(bad_sisters) > 0 or len(bad_subsisters) > 0:
			#This loop will keep moving towards the root of the tree until it finds a node that
			#has leaves from a cell other than the one for which we are looking for sisters
			parent_node = leaf; seen_taxa = {coconts[leaf.name[:10]]}; sisters = []
			while len(seen_taxa) == 1:
				parent_node = parent_node.up

				for l2 in parent_node:
					seen_taxa.add(coconts[l2.name[:10]])
					if coconts[l2.name[:10]] != coconts[leaf.name[:10]]:
						sisters.append(l2.name[:10])

			#Create a record of the subsister sequences
			sub_sisters = []
			if args.subsister_rules != None:
				new_parent_node = parent_node.up

				for l2 in new_parent_node:
					if l2.name not in [l.name for l in parent_node]:
						sub_sisters.append(l2.name)

			#Create a record of the sister sequences
			sisters = list(dict.fromkeys(sisters))

			#Getting list of removable sequences by sister relationships
			sisters_removable = []; bls = []
			for contam in bad_sisters:
				for sister in sisters:
					if sister.startswith(contam) and sister not in sisters_removable:
						sisters_removable.append(sister)
						bls.append(bad_sisters[contam])

			#Getting list of removable sequences by sub-sister relationships
			subsisters_removable = []
			for contam in bad_subsisters:
				for sub_sister in sub_sisters:
					if sub_sister.startswith(contam) and sub_sister not in subsisters_removable:
						subsisters_removable.append(sub_sister)

			if len(bls) > 0:
				bl_rule_min = min(bls)
			else:
				bl_rule_min = 0

			if len(sisters_removable) == len(sisters) and leaf.dist <= bl_rule_min*mean_bl and len(sisters_removable) > 0:
				seqs2remove.append(leaf.name)
			elif len(subsisters_removable) == len(sub_sisters) and len(subsisters_removable) > 0:
				seqs2remove.append(leaf.name)

	return [leaf.name for leaf in tree if leaf.name not in seqs2remove]

def write_new_preguidance(params, seqs2keep, seqs_per_og, tree_file):
	"""Write a new unaligned sequence file excluding removed sequences.

	Match the tree file to its corresponding sequence file and write the sequences that should
	be kept to Output/Pre-Guidance. Also exclude taxa listed in params.cl_exclude_taxa.

	Args:
		params (argparse.Namespace): Command-line arguments.
		seqs2keep (list[str]): Names of sequences to keep.
		seqs_per_og (dict): Sequences grouped by input file, in the form {file name: {sequence name: sequence}}.
		tree_file (str): Name of the tree file.

	Returns:
		tuple: The matching sequence file name and a list of removed sequence names.
		Returns (None, []) if no unique matching sequence file is found.
	"""

	if params.cl_exclude_taxa != None:
		try:
			exclude_taxa = list(dict.fromkeys([line.strip() for line in open(params.cl_exclude_taxa)]))
		except (FileNotFoundError, TypeError) as e:
			print('\nERROR: Unable to read the file listing taxa to exclude in the first iteration of the contamination loop (--cl_exclude_taxa). Please make sure that the path is correct and that the file is formatted correctly.\n\n' + str(e) + '\n')
			exit()
	else:
		exclude_taxa = []

	prefix = tree_file.split('.')[0]
	seq_file = [file for file in seqs_per_og if file.startswith(prefix)]

	if len(seq_file) == 0:
		print('\nNo sequence file found for tree file ' + tree_file + '. Skipping this gene family.\n')
		return None, []
	elif len(seq_file) > 1:
		print('\nMore than one sequence file found matching the tree file ' + tree_file + '. Please make your file names more unique: there should be one sequence file for every tree file, with a matching unique prefix (everything before the first "."). Skipping this gene family.\n')
		return None, []
	elif len(seq_file) == 1:
		with open(params.output + '/Output/Pre-Guidance/' + seq_file[0], 'w') as o:
			for rec in seqs_per_og[seq_file[0]]:
				if rec in seqs2keep and rec[:10] not in exclude_taxa and rec[:2] not in exclude_taxa and rec[:5] not in exclude_taxa:
					o.write('>' + rec + '\n' + seqs_per_og[seq_file[0]][rec] + '\n\n')

		seqs_removed_from_og = [seq for seq in seqs_per_og[seq_file[0]] if seq not in seqs2keep]

		return seq_file[0], seqs_removed_from_og


def cl_mafft(params):
	"""Run MAFFT on all FASTA files in the Pre-Guidance directory and save the untrimmed alignments to the NotGapTrimmed directory. Then, use trimAl to remove gaps based on the specified gap threshold and save the trimmed alignments to the Guidance directory.

	Args:
		params (argparse.Namespace): Command-line arguments containing the output directory and trimAl gap threshold (trimal_cutoff).
	"""
	for file in os.listdir(params.output + '/Output/Pre-Guidance'):
		if file.split('.')[-1] in ('fasta', 'fas', 'faa'):
			os.system('mafft ' + params.output + '/Output/Pre-Guidance/' + file + ' > ' + params.output + '/Output/NotGapTrimmed/' + file)

			os.system('Scripts/trimal-trimAl/source/trimal -in ' + params.output + '/Output/NotGapTrimmed/' + file + ' -out ' + params.output + '/Output/Guidance/' + file.split('.')[0] + '.70gapTrimmed.fasta' + ' -gapthreshold ' + str(params.trimal_cutoff) + ' -fasta')

def cl_fasttree(params):
	"""Build phylogenetic trees using FastTree.

	Run FastTree on all FASTA files in the Guidance directory and save the resulting trees to the Trees directory,
	the Output/Trees under the base directory specified by params.output.

	Args:
		params (argparse.Namespace): Command-line arguments containing the output directory path.
	"""
	for file in os.listdir(params.output + '/Output/Guidance'):
		if file.split('.')[-1] in ('fasta', 'fas', 'faa'):
			os.system('FastTree ' + params.output + '/Output/Guidance/' + file + ' > ' + params.output + '/Output/Trees/' + file.split('.')[0] + '.FastTree.tre')


def cl_iqtree(params):
	"""Build phylogenetic trees using IQ-TREE.

	Run Iqtree on all FASTA files in the Guidance directory and save the resulting trees to the Trees directory,
	the Output/Trees under the base directory specified by params.output.

	Args:
		params (argparse.Namespace): Command-line arguments containing the output directory path.
	"""
	for file in os.listdir(params.output + '/Output/Guidance'):
		if file.split('.')[-1] in ('fasta', 'fas', 'faa'):
			if not os.path.isdir(params.output + '/Output/Intermediate/IQTree'):
				os.mkdir(params.output + '/Output/Intermediate/IQTree')
			tax_iqtree_outdir = params.output + '/Output/Intermediate/IQTree/' + file.split('.')[0].split('_preguidance')[0]
			os.mkdir(tax_iqtree_outdir)
			os.system('iqtree2 -s ' + params.output + '/Output/Guidance/' + file + ' -m LG+G -T 10 --prefix ' +  tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree')
			#Copy over the final output
			if os.path.isfile(tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.treefile'):
				os.system('cp ' + tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.treefile ' + params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.tree')
	os.system('rm -r ' + params.output + '/Output/Intermediate/IQTree/*')


def cl_iqtree_fast(params):
	"""Build phylogenetic trees using IQ-TREE in fast mode.

	Run Iqtree_fast on all FASTA files in the Guidance directory and save the resulting trees to the Trees directory,
	the Output/Trees under the base directory specified by params.output.

	Args:
		params (argparse.Namespace): Command-line arguments containing the output directory path.
	"""
	for file in os.listdir(params.output + '/Output/Guidance'):
		if file.split('.')[-1] in ('fasta', 'fas', 'faa'):
			if not os.path.isdir(params.output + '/Output/Intermediate/IQTree'):
				os.mkdir(params.output + '/Output/Intermediate/IQTree')
			tax_iqtree_outdir = params.output + '/Output/Intermediate/IQTree/' + file.split('.')[0].split('_preguidance')[0]
			os.mkdir(tax_iqtree_outdir)
			os.system('iqtree2 -s ' + params.output + '/Output/Guidance/' + file + ' -m LG+G -T 10 --fast --prefix ' +  tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree')
			#Copy over the final output
			if os.path.isfile(tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.treefile'):
				os.system('cp ' + tax_iqtree_outdir + '/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.treefile ' + params.output + '/Output/Trees/' + file.split('.')[0].split('_preguidance')[0] + '.IQTree.tree')
	os.system('rm -r ' + params.output + '/Output/Intermediate/IQTree/*')


def run(params):
	"""Run the iterative contamination removal pipeline.

	Repeat the process params.nloops times:
	1) Read sequence files from Output/Pre-Guidance when starting from raw data, or from params.data otherwise.
	2) For each tree not yet completed, identify sequences to keep using clade-based (params.contamination_loop == 'clade') or sister-based
		(params.contamination_loop == 'seq') removal, then write a new unaligned sequence file. Mark the tree as completed when the returned list of removed sequences is empty.
	3) Record removed sequences and their iteration numbers in params.output/Output/SequencesRemoved_ContaminationLoop.txt.
	4) Archive the Trees, Guidance, and NotGapTrimmed directories, appending the iteration number to each directory name.
	5) Re-align sequences using MAFFT or Guidance according to params.cl_alignment_method, then rebuild trees according to params.cl_tree_method.


    Args:
        params (argparse.Namespace): Command-line arguments specifying the number of iterations, contamination removal method, input and
        	output directories, alignment method, tree-building method, and optional contamination rules.

    Notes:
        The function modifies params.start, params.end, and params.tree_method during execution to configure the alignment and tree-building steps.

		The function exits if the specified sub-sister rules file cannot be found.
	"""
	seqs_removed = []
	completed_ogs = []

	with open('SequencesRemoved_ContaminationLoop.txt', 'w') as o:
		o.write('Sequence\tLoopRemoved\n')

	#For each iteration
	for loop in range(params.nloops):
		seqs_removed_loop = []

		#Finding input files
		if params.start == 'raw':
			seqs_per_og = { file : { rec.id : str(rec.seq) for rec in SeqIO.parse(params.output + '/Output/Pre-Guidance/' + file, 'fasta') } for file in os.listdir(params.output + '/Output/Pre-Guidance') if file.split('.')[-1] in ('fasta', 'fas', 'faa') }
		elif params.start in ('unaligned', 'aligned', 'trees'):
			seqs_per_og = { file : { rec.id : str(rec.seq).replace('-', '') for rec in SeqIO.parse(params.data + '/' + file, 'fasta') } for file in os.listdir(params.data) if file.split('.')[-1] in ('fasta', 'fas', 'faa') }

			if loop == 0:
				for file in os.listdir(params.data):
					if file.split('.')[-1] in ('tre', 'tree', 'treefile'):
						os.system('cp ' + params.data + '/' + file + ' ' + params.output + '/Output/Trees')
		if loop > 0 or params.start == 'raw':
			os.system('mv ' + params.output + '/Output/Pre-Guidance ' + params.output + '/Output/Pre-Guidance_' + str(loop))

			os.mkdir(params.output + '/Output/Pre-Guidance')

		#Wrapper for running clade-based contamination removal on all trees
		if params.contamination_loop == 'clade':
			for tree_file in os.listdir(params.output + '/Output/Trees'):
				if tree_file.split('.')[-1] in ('tre', 'tree', 'treefile', 'nex') and tree_file not in completed_ogs:
					seqs2keep = get_subtrees(params, params.output + '/Output/Trees/' + tree_file)

					seq_file, seqs_removed_from_og = write_new_preguidance(params, seqs2keep, seqs_per_og, tree_file)

					if len(seqs_removed_from_og) == 0:
						completed_ogs.append(tree_file)
					else:
						seqs_removed_loop += [seq for seq in seqs_per_og[seq_file] if seq not in seqs2keep and seq not in seqs_removed]
		#Wrapper for running sisters-based contamination removal on all trees
		elif params.contamination_loop == 'seq':

			sister_contam_per_tax = { }
			if params.sister_rules != None:
				for line in open(params.sister_rules):
					if line.strip().split('\t')[0] not in sister_contam_per_tax:
						sister_contam_per_tax.update({ line.strip().split('\t')[0] : [] })

					try:
						sister_contam_per_tax[line.strip().split('\t')[0]].append((line.strip().split('\t')[1], float(line.strip().split('\t')[2])))
					except ValueError:
						sister_contam_per_tax[line.strip().split('\t')[0]].append((line.strip().split('\t')[1], float('inf')))
					except IndexError:
						if line.strip() != '':
							print('\nWarning: the line "' + line.strip() + '" in the sister rules file could not be processed\n')

			subsister_contam_per_tax = { }
			if params.subsister_rules != None:
				if os.path.isfile(params.subsister_rules):
					for line in open(params.subsister_rules):
						if line.strip().split('\t')[0] not in subsister_contam_per_tax and len(line.strip().split('\t')) == 2:
							subsister_contam_per_tax.update({ line.strip().split('\t')[0] : [] })

						subsister_contam_per_tax[line.strip().split('\t')[0]].append(line.strip().split('\t')[1])
				else:
					print('\nERROR: It looks like you tried to input a sub-sister rules file to the contamination loop, but the file could not be found.\n')
					exit()

			for tree_file in os.listdir(params.output + '/Output/Trees'):
				if tree_file.split('.')[-1] in ('tre', 'tree', 'treefile', 'nex') and tree_file not in completed_ogs:
					seqs2keep = get_sisters(params, params.output + '/Output/Trees/' + tree_file, sister_contam_per_tax, subsister_contam_per_tax)

					seq_file, seqs_removed_from_og = write_new_preguidance(params, seqs2keep, seqs_per_og, tree_file)

					if len(seqs_removed_from_og) == 0:
						completed_ogs.append(tree_file)
					else:
						seqs_removed_loop += [seq for seq in seqs_per_og[seq_file] if seq not in seqs2keep and seq not in seqs_removed]

		#Keeping record of removed sequences
		seqs_removed += seqs_removed_loop
		with open(params.output + '/Output/SequencesRemoved_ContaminationLoop.txt', 'a') as o:
			for seq in seqs_removed_loop:
				o.write(seq + '\t' + str(loop) + '\n')

		#Writing output files with sequences removed, with the iteration labeled
		os.system('mv ' + params.output + '/Output/Trees ' + params.output + '/Output/Trees_' + str(loop))
		os.mkdir(params.output + '/Output/Trees')

		os.system('mv ' + params.output + '/Output/Guidance ' + params.output + '/Output/Guidance_' + str(loop))
		os.mkdir(params.output + '/Output/Guidance')

		os.system('mv ' + params.output + '/Output/NotGapTrimmed ' + params.output + '/Output/NotGapTrimmed_' + str(loop))
		os.mkdir(params.output + '/Output/NotGapTrimmed')

		params.start = 'unaligned'
		params.end = 'trees'
		params.tree_method = params.cl_tree_method

		#Re-aligning and building trees without contaminant sequences... then ready for next iteration.
		if params.cl_alignment_method == 'mafft_only':
			cl_mafft(params)
		else:
			guidance.run(params)

		if params.cl_tree_method == 'fasttree':
			cl_fasttree(params)
		elif params.cl_tree_method == 'iqtree':
			cl_iqtree(params)
			os.system('rm -r ' + params.output + '/Output/Intermediate/IQTree/*')
		elif params.cl_tree_method == 'iqtree_fast':
			cl_iqtree_fast(params)
			os.system('rm -r ' + params.output + '/Output/Intermediate/IQTree/*')
		elif params.cl_tree_method == 'raxml':
			os.system('rm -r ' + params.output + '/Output/Intermediate/RAxML/*')