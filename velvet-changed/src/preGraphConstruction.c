/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "globals.h"
#include "preGraph.h"
#include "recycleBin.h"
#include "roadMap.h"
#include "readSet.h"
#include "concatenatedPreGraph.h"
#include "utility.h"
#include "kmer.h"
#include "tightString.h"
#include "binarySequences.h"
#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

#ifdef _OPENMP

Coordinate *annotationOffset = NULL;

static omp_lock_t *nodeLocks = NULL;

static void createNodeLocks(PreGraph *preGraph)
{
	IDnum nbNodes;
	IDnum nodeIndex;

	nbNodes = preNodeCount_pg(preGraph) + 1;
	if (nodeLocks)
		free (nodeLocks);
	nodeLocks = mallocOrExit(nbNodes, omp_lock_t);

	#pragma omp parallel for
	for (nodeIndex = 0; nodeIndex < nbNodes; nodeIndex++)
		omp_init_lock(nodeLocks + nodeIndex);
}

static void lockNode(IDnum preNodeID)
{
	omp_set_lock(nodeLocks + preNodeID);
}

static void unLockNode(IDnum preNodeID)
{
	omp_unset_lock(nodeLocks + preNodeID);
}

static void lockTwoNodes(IDnum preNodeID, IDnum preNode2ID)
{
	if (preNodeID < 0)
		preNodeID = -preNodeID;
	if (preNode2ID < 0)
		preNode2ID = -preNode2ID;

	/* Lock lowest ID first to avoid deadlocks */
	if (preNodeID == preNode2ID)
		omp_set_lock (nodeLocks + preNodeID);
	else if (preNodeID < preNode2ID)
	{
		omp_set_lock (nodeLocks + preNodeID);
		omp_set_lock (nodeLocks + preNode2ID);
	}
	else
	{
		omp_set_lock (nodeLocks + preNode2ID);
		omp_set_lock (nodeLocks + preNodeID);
	}
}

static void unLockTwoNodes(IDnum preNodeID, IDnum preNode2ID)
{
	if (preNodeID < 0)
		preNodeID = -preNodeID;
	if (preNode2ID < 0)
		preNode2ID = -preNode2ID;

	omp_unset_lock (nodeLocks + preNodeID);
	if (preNodeID != preNode2ID)
		omp_unset_lock (nodeLocks + preNode2ID);
}
#endif

// Internal structure used to mark the ends of an Annotation
struct insertionMarker_st {
	Annotation *annot;
	boolean isStart;
}  ATTRIBUTE_PACKED;

Coordinate getInsertionMarkerPosition(InsertionMarker * marker)
{
	if (marker->isStart)
		return getStart(marker->annot);
	else
		return getFinish(marker->annot);
}

int compareInsertionMarkers(const void *A, const void *B)
{
	Coordinate Apos =
	    getInsertionMarkerPosition((InsertionMarker *) A);
	Coordinate Bpos =
	    getInsertionMarkerPosition((InsertionMarker *) B);

	if (Apos < Bpos)
		return -1;
	else if (Apos == Bpos)
		return 0;
	else
		return 1;
}

// Applies mergeSort to each insertion marker list (in order of position)
static void
orderInsertionMarkers(InsertionMarker ** insMarkers,
		      IDnum * markerCounters, RoadMapArray * rdmaps)
{
	IDnum sequenceIndex;
	IDnum sequenceCounter = rdmaps->length;

	velvetLog("Ordering insertion markers\n");
#ifdef _OPENMP
	#pragma omp parallel for
#endif
	for (sequenceIndex = 1; sequenceIndex <= sequenceCounter;
	     sequenceIndex++) {
		qsort(insMarkers[sequenceIndex],
		      markerCounters[sequenceIndex],
		      sizeof(InsertionMarker), compareInsertionMarkers);
	}
}

// Creates insertion marker lists 
static void
setInsertionMarkers(RoadMapArray * rdmaps,
		    IDnum * markerCounters,
		    InsertionMarker ** veryLastMarker,
		    InsertionMarker ** insertionMarkers)
{
	IDnum sequenceCounter = rdmaps->length;
	IDnum sequenceIndex, sequenceIndex2;
	Coordinate totalCount = 0;
	RoadMap *rdmap;
	Annotation *annot = rdmaps->annotations;
	InsertionMarker *nextMarker, *newMarker;
	IDnum annotIndex, lastAnnotIndex;
	InsertionMarker **insMarkers =
	    callocOrExit(rdmaps->length + 1, InsertionMarker *);
	// Counting insertion markers
	for (sequenceIndex = 1; sequenceIndex < sequenceCounter + 1;
	     sequenceIndex++) {
		//velvetLog("Going through sequence %d\n", sequenceIndex);
		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		lastAnnotIndex = getAnnotationCount(rdmap);

		// Set insertion markers in previous sequences :

		for (annotIndex = 0; annotIndex < lastAnnotIndex;
		     annotIndex++) {
			if (getAnnotSequenceID(annot) > 0) {
				markerCounters[getAnnotSequenceID(annot)]
				    += 2;
			} else {
				markerCounters[-getAnnotSequenceID(annot)]
				    += 2;
			}
			totalCount += 2;
			annot = getNextAnnotation(annot);
		}
	}

	// Allocating space
	*insertionMarkers = callocOrExit(totalCount, InsertionMarker);
	*veryLastMarker = *insertionMarkers + totalCount;

	// Pointing each node to its space      
	nextMarker = *insertionMarkers;
	for (sequenceIndex = 1; sequenceIndex < sequenceCounter + 1;
	     sequenceIndex++) {
		insMarkers[sequenceIndex] = nextMarker;
		nextMarker = nextMarker + markerCounters[sequenceIndex];
		markerCounters[sequenceIndex] = 0;
	}

	// Filling up space with data
	annot = rdmaps->annotations;
	for (sequenceIndex = 1; sequenceIndex < sequenceCounter + 1;
	     sequenceIndex++) {
		//velvetLog("Going through sequence %d\n", sequenceIndex);
		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		lastAnnotIndex = getAnnotationCount(rdmap);

		// Set insertion markers in previous sequences :

		for (annotIndex = 0; annotIndex < lastAnnotIndex;
		     annotIndex++) {
			sequenceIndex2 = getAnnotSequenceID(annot);
			if (sequenceIndex2 > 0) {
				newMarker =
				    insMarkers[sequenceIndex2] +
				    (markerCounters[sequenceIndex2])++;
				newMarker->annot = annot;
				newMarker->isStart = true;

				newMarker =
				    insMarkers[sequenceIndex2] +
				    (markerCounters[sequenceIndex2])++;
				newMarker->annot = annot;
				newMarker->isStart = false;
			} else {
				incrementAnnotationCoordinates(annot);

				newMarker =
				    insMarkers[-sequenceIndex2] +
				    (markerCounters[-sequenceIndex2])++;
				newMarker->annot = annot;
				newMarker->isStart = true;

				newMarker =
				    insMarkers[-sequenceIndex2] +
				    (markerCounters[-sequenceIndex2])++;
				newMarker->annot = annot;
				newMarker->isStart = false;
			}
			annot = getNextAnnotation(annot);
		}
	}

	orderInsertionMarkers(insMarkers, markerCounters, rdmaps);
	free(insMarkers);
}

// Counts how many preNodes are to be created to allocate appropriate memory
static void
countPreNodes(RoadMapArray * rdmaps, PreGraph * preGraph,
	      IDnum * markerCounters, InsertionMarker * insertionMarkers,
	      InsertionMarker * veryLastMarker)
{
	Annotation *annot = rdmaps->annotations;
	InsertionMarker *currentMarker = insertionMarkers;
	IDnum markerIndex, lastMarkerIndex;
	IDnum sequenceIndex;
	Coordinate currentPosition, nextStop;
	IDnum preNodeCounter = 0;
	RoadMap *rdmap;
	IDnum annotIndex, lastAnnotIndex;


	// Now that we have read all of the annotations, we go on to create the preNodes and tie them up
	for (sequenceIndex = 1;
	     sequenceIndex <= sequenceCount_pg(preGraph);
	     sequenceIndex++) {
		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		annotIndex = 0;
		lastAnnotIndex = getAnnotationCount(rdmap);
		markerIndex = 0;
		lastMarkerIndex = markerCounters[sequenceIndex];
		currentPosition = 0;

		while (annotIndex < lastAnnotIndex) {
			if (markerIndex == lastMarkerIndex
			    || getPosition(annot) <=
			    getInsertionMarkerPosition(currentMarker))
				nextStop = getPosition(annot);
			else
				nextStop =
				    getInsertionMarkerPosition
				    (currentMarker);

			if (currentPosition != nextStop) {
				preNodeCounter++;
				currentPosition = nextStop;
			}

			while (markerIndex < lastMarkerIndex
			       && getInsertionMarkerPosition(currentMarker)
			       == currentPosition) {
				currentMarker++;
				markerIndex++;
			}

			while (annotIndex < lastAnnotIndex
			       && getPosition(annot) == currentPosition) {
				annot = getNextAnnotation(annot);
				annotIndex++;
			}

		}

		while (markerIndex < lastMarkerIndex) {
			if (currentPosition ==
			    getInsertionMarkerPosition(currentMarker)) {
				currentMarker++;
				markerIndex++;
			} else {
				preNodeCounter++;
				currentPosition =
				    getInsertionMarkerPosition
				    (currentMarker);
			}
		}
	}

	allocatePreNodeSpace_pg(preGraph, preNodeCounter);
}

// Counts how many preNodes are to be created to allocate appropriate memory
static void
countPreNodesOMP(RoadMapArray * rdmaps, PreGraph * preGraph,
	      IDnum * markerCounters, InsertionMarker * insertionMarkers,
	      InsertionMarker * veryLastMarker, IDnum *chains)
{
	Annotation *annot = rdmaps->annotations;
	InsertionMarker *currentMarker = insertionMarkers;
	IDnum markerIndex, lastMarkerIndex;
	IDnum sequenceIndex;
	Coordinate currentPosition, nextStop;
	IDnum preNodeCounter = 0;
	RoadMap *rdmap;
	IDnum annotIndex, lastAnnotIndex;

	// IDnum *markerOffsets = malloc((sequenceCount_pg(preGraph) + 1) * sizeof(IDnum));
	// markerOffsets[0] = 0;
	// for (IDnum i = 1; i <= sequenceCount_pg(preGraph); i++) {
	// 	markerOffsets[i] = markerOffsets[i - 1] + markerCounters[i];
	// }
	
	// Get starting annotation for all sequences with a prefix sum
	Annotation **startAnnot = malloc((sequenceCount_pg(preGraph) + 1) * sizeof(Annotation *));
	startAnnot[0] = annot;
	for (IDnum i = 1; i <= sequenceCount_pg(preGraph); i++) {
		startAnnot[i] = startAnnot[i - 1] + getAnnotationCount(getRoadMapInArray(rdmaps, i - 1));
	}

	// Get starting marker for all sequences with a prefix sum
	InsertionMarker **startMarkers = malloc((sequenceCount_pg(preGraph) + 1) * sizeof(InsertionMarker *));
	startMarkers[0] = currentMarker;
	for (IDnum i = 1; i <= sequenceCount_pg(preGraph); i++) {
		startMarkers[i] = startMarkers[i - 1] + markerCounters[i];
		// startMarkers[i] = startMarkers[i - 1] + markerCounters[i - 1];

	}

	IDnum *preNodeCounts = malloc((sequenceCount_pg(preGraph) + 1) * sizeof(IDnum));  // N = number of sequences

	// Now that we have read all of the annotations, we go on to create the preNodes and tie them up
	velvetLog("Number of sequences = %li\n",
		(long) sequenceCount_pg(preGraph));
 	#pragma omp parallel for private(rdmap, annotIndex, lastAnnotIndex, annot, markerIndex, lastMarkerIndex, currentPosition, currentMarker) reduction(+:preNodeCounter)
	for (sequenceIndex = 1;
	     sequenceIndex <= sequenceCount_pg(preGraph);
	     sequenceIndex++) {
		// InsertionMarker *markerStart = &insertionMarkers[markerOffsets[sequenceIndex - 1]];

		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		annot = startAnnot[sequenceIndex - 1];
		annotIndex = 0;
		lastAnnotIndex = getAnnotationCount(rdmap);
		markerIndex = 0;
		lastMarkerIndex = markerCounters[sequenceIndex];
		currentPosition = 0;
		currentMarker = startMarkers[sequenceIndex - 1];
		IDnum preNodeCount = 0;


		while (annotIndex < lastAnnotIndex) {
			Coordinate nextStop;
			if (markerIndex == lastMarkerIndex
			    || getPosition(annot) <=
			    getInsertionMarkerPosition(currentMarker))
				nextStop = getPosition(annot);
			else
				nextStop = getInsertionMarkerPosition(currentMarker);

			if (currentPosition != nextStop) {
				preNodeCount++;
				preNodeCounter++;
				currentPosition = nextStop;
			}

			while (markerIndex < lastMarkerIndex
			       && getInsertionMarkerPosition(currentMarker)
			       == currentPosition) {
				currentMarker++;
				markerIndex++;
			}

			while (annotIndex < lastAnnotIndex
			       && getPosition(annot) == currentPosition) {
				annot = getNextAnnotation(annot); // annot += 1
				annotIndex++;
			}

		}

		while (markerIndex < lastMarkerIndex) {
			if (currentPosition ==
			    getInsertionMarkerPosition(currentMarker)) {
				currentMarker++;
				markerIndex++;
			} else {
				preNodeCount++;
				preNodeCounter++;
				currentPosition = getInsertionMarkerPosition(currentMarker);
			}
		}

		preNodeCounts[sequenceIndex] = preNodeCount;
	}

	chains[1] = 1;  // or 0 if you're 0-based
	for (int i = 2; i <= sequenceCount_pg(preGraph); ++i)
    	chains[i] = chains[i - 1] + preNodeCounts[i - 1];

	velvetLog("preNodeCounter = %li, chains[N] + preNodeCounts[N] = %li\n",
		(long) preNodeCounter,
		(long) chains[sequenceCount_pg(preGraph)] +
		preNodeCounts[sequenceCount_pg(preGraph)]);
	// free(markerOffsets);
	free(preNodeCounts);
	free(startMarkers);
	free(startAnnot);
	allocatePreNodeSpace_pg(preGraph, preNodeCounter);
}

static void convertInsertionMarkers(InsertionMarker * insertionMarkers,
				    InsertionMarker * veryLastMarker,
				    IDnum * chains)
{
	InsertionMarker *marker;
	Annotation *annot;

	for (marker = insertionMarkers; marker != veryLastMarker; marker++) {
		annot = marker->annot;

		if (getAnnotSequenceID(annot) > 0) {
			if (marker->isStart) {
				if (getStartID(annot) == 0)
					setStartID(annot,
						   chains
						   [getAnnotSequenceID
						    (annot)]);
				else
					setStartID(annot,
						   getStartID(annot) + 1);
			}
		} else {
			if (marker->isStart)
				setStartID(annot, -getStartID(annot));
			else {
				if (getFinishID(annot) == 0)
					setFinishID(annot,
						    -chains
						    [-getAnnotSequenceID
						     (annot)]);
				else
					setFinishID(annot,
						    -getFinishID(annot) -
						    1);
			}
		}
	}

	free(insertionMarkers);
}

static void convertMarker(InsertionMarker * marker, IDnum nodeID)
{
	if (marker->isStart)
		setStartID(marker->annot, nodeID);
	else
		setFinishID(marker->annot, nodeID);
}


// Creates the preNode using insertion marker and annotation lists for each sequence
static void
// Creates the preNode using insertion marker and annotation lists for each sequence
createPreNodes(RoadMapArray * rdmaps, PreGraph * preGraph,
	       IDnum * markerCounters, InsertionMarker * insertionMarkers,
	       InsertionMarker * veryLastMarker, IDnum * chains,
	       SequencesReader *seqReadInfo, int WORDLENGTH)
{
	char *sequenceFilename = seqReadInfo->m_seqFilename;
	Annotation *annot = rdmaps->annotations;
	IDnum latestPreNodeID;
	InsertionMarker *currentMarker = insertionMarkers;
	IDnum sequenceIndex;
	Coordinate currentPosition, nextStop;
	IDnum preNodeCounter = 1;
	FILE *file = NULL;
	char line[50000];
	int lineLength = 50000;
	Coordinate readIndex;
	boolean tooShort;
	Kmer initialKmer;
	char c;
	RoadMap *rdmap;
	IDnum annotIndex, lastAnnotIndex;
	IDnum markerIndex, lastMarkerIndex;

	if (!seqReadInfo->m_bIsBinary) {
		file = fopen(sequenceFilename, "r");
	if (file == NULL) 
		exitErrorf(EXIT_FAILURE, true, "Could not read %s", sequenceFilename);
	// Reading sequence descriptor in first line
	if (sequenceCount_pg(preGraph) > 0 && !fgets(line, lineLength, file))
		exitErrorf(EXIT_FAILURE, true, "%s incomplete.", sequenceFilename);
		seqReadInfo->m_pFile = file;
	}

	// Now that we have read all of the annotations, we go on to create the preNodes and tie them up
	for (sequenceIndex = 1;
	     sequenceIndex <= sequenceCount_pg(preGraph);
	     sequenceIndex++) {
		if (sequenceIndex % 1000000 == 0)
			velvetLog("Sequence %li / %li\n", (long) sequenceIndex,
			       (long) sequenceCount_pg(preGraph));

		if (!seqReadInfo->m_bIsBinary) {
		while (line[0] != '>')
			if (!fgets(line, lineLength, file))
				exitErrorf(EXIT_FAILURE, true, "%s incomplete.", sequenceFilename);
		}

		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		annotIndex = 0;
		lastAnnotIndex = getAnnotationCount(rdmap);
		markerIndex = 0;
		lastMarkerIndex = markerCounters[sequenceIndex];
		currentPosition = 0;

		// Reading first (k-1) nucleotides
		tooShort = false;
		clearKmer(&initialKmer);
		//velvetLog("Initial kmer: ");
		TightString *tString = NULL;
		char *strString = NULL;
		if (seqReadInfo->m_bIsBinary) {
			tString = getTightStringInArray(seqReadInfo->m_sequences->tSequences, sequenceIndex - 1);
			strString = readTightString(tString);
		}
		for (readIndex = 0; readIndex < WORDLENGTH - 1;
		     readIndex++) {
			if (seqReadInfo->m_bIsBinary) {
				if (readIndex >= tString->length) {
					tooShort = true;
					break;
				}

				c = strString[readIndex];
			} else {
			c = getc(file);
			while (c == '\n' || c == '\r') 
				c = getc(file);
	
			if (c == '>' || c == 'M' || c == EOF) {
				ungetc(c, file);
				tooShort = true;
				break;
			}
			}
			switch (c) {
			case 'A':
			case 'N':
				pushNucleotide(&initialKmer, ADENINE);
				break;
			case 'C':
				pushNucleotide(&initialKmer, CYTOSINE);
				break;
			case 'G':
				pushNucleotide(&initialKmer, GUANINE);
				break;
			case 'T':
				pushNucleotide(&initialKmer, THYMINE);
				break;
			default:
				velvetLog
				    ("Irregular sequence file: are you sure your Sequence and Roadmap file come from the same source?\n");
				fflush(stdout);
				abort();
			}
		}

		if (tooShort) {
			//velvetLog("Skipping short read.. %d\n", sequenceIndex);
			chains[sequenceIndex] = preNodeCounter;
			if (seqReadInfo->m_bIsBinary) {
				free(strString);
			} else {
			if (!fgets(line, lineLength, file) && sequenceIndex < sequenceCount_pg(preGraph))
				exitErrorf(EXIT_FAILURE, true, "%s incomplete.", sequenceFilename);
			}
			continue;
		}

		char *currString = NULL;
		if (seqReadInfo->m_bIsBinary) {
			currString = &strString[readIndex];
			seqReadInfo->m_ppCurrString = &currString;
		}
		latestPreNodeID = 0;

		while (annotIndex < lastAnnotIndex) {
			if (markerIndex == lastMarkerIndex
			    || getPosition(annot) <=
			    getInsertionMarkerPosition(currentMarker))
				nextStop = getPosition(annot);
			else {
				nextStop =
				    getInsertionMarkerPosition
				    (currentMarker);
			}

			if (currentPosition != nextStop) {
				if (seqReadInfo->m_bIsBinary) {
					if (readIndex >= tString->length) {
						velvetLog("readIndex %ld beyond string len %ld\n", (uint64_t) readIndex, (uint64_t) tString->length);
						exit(1);
					}
				}
				//if (sequenceIndex == 481)
				//	velvetLog("Adding pre nodes from %lli to %lli\n", (long long) currentPosition, (long long) nextStop);
				addPreNodeToPreGraph_pg(preGraph,
							currentPosition,
							nextStop,
							seqReadInfo,
							&initialKmer,
							preNodeCounter);
				if (latestPreNodeID == 0) {
					chains[sequenceIndex] =
					    preNodeCounter;
				}
				latestPreNodeID = preNodeCounter++;
				currentPosition = nextStop;
			}

			while (markerIndex < lastMarkerIndex
			       && getInsertionMarkerPosition(currentMarker)
			       == nextStop) {
				convertMarker(currentMarker,
					      latestPreNodeID);
				currentMarker++;
				markerIndex++;
			}

			while (annotIndex < lastAnnotIndex
			       && getPosition(annot) == nextStop) {
				for (readIndex = 0;
				     readIndex <
				     getAnnotationLength(annot);
				     readIndex++) {
					if (seqReadInfo->m_bIsBinary) {
						c = *currString;
						currString += 1;   // increment the pointer
					} else {
					c = getc(file);
					while (!isalpha(c))
						c = getc(file);
					}

					//if (sequenceIndex == 481)
					//	velvetLog("(%c)", c);
					switch (c) {
					case 'A':
					case 'N':
						pushNucleotide(&initialKmer, ADENINE);
						break;
					case 'C':
						pushNucleotide(&initialKmer, CYTOSINE);
						break;
					case 'G':
						pushNucleotide(&initialKmer, GUANINE);
						break;
					case 'T':
						pushNucleotide(&initialKmer, THYMINE);
						break;
					default:
						velvetLog
						    ("Irregular sequence file: are you sure your Sequence and Roadmap file come from the same source?\n");
						fflush(stdout);
#ifdef DEBUG 
						abort();
#endif 
						exit(1);
					}
				}

				annot = getNextAnnotation(annot);
				annotIndex++;
			}

		}

		while (markerIndex < lastMarkerIndex) {
			if (currentPosition ==
			    getInsertionMarkerPosition(currentMarker)) {
				convertMarker(currentMarker,
					      latestPreNodeID);
				currentMarker++;
				markerIndex++;
			} else {
				nextStop =
				    getInsertionMarkerPosition
				    (currentMarker);
				//if (sequenceIndex == 481)
				//	velvetLog("Adding pre nodes from %lli to %lli\n", (long long) currentPosition, (long long) nextStop);
				addPreNodeToPreGraph_pg(preGraph,
							currentPosition,
							nextStop, seqReadInfo,
							&initialKmer,
							preNodeCounter);
				if (latestPreNodeID == 0)
					chains[sequenceIndex] =
					    preNodeCounter;
				latestPreNodeID = preNodeCounter++;
				currentPosition =
				    getInsertionMarkerPosition
				    (currentMarker);
			}
		}
		if (seqReadInfo->m_bIsBinary) {
			free(strString);
		} else {
		// End of sequence
		if (!fgets(line, lineLength, file) && sequenceIndex < sequenceCount_pg(preGraph))
			exitErrorf(EXIT_FAILURE, true, "%s incomplete.", sequenceFilename);
		//velvetLog(" \n");
		}

		if (latestPreNodeID == 0)
			chains[sequenceIndex] = preNodeCounter;
	}

	velvetLog("preNodeCounter = %li, actual number of preNodes = %li\n", (long) preNodeCounter, (long) preNodeCount_pg(preGraph));
	free(markerCounters);
	if (!seqReadInfo->m_bIsBinary) {
	fclose(file);
	}

}

// Creates the preNode using insertion marker and annotation lists for each sequence
static void
// Creates the preNode using insertion marker and annotation lists for each sequence
createPreNodesOMP(RoadMapArray * rdmaps, PreGraph * preGraph,
	       IDnum * markerCounters, InsertionMarker * insertionMarkers,
	       InsertionMarker * veryLastMarker, IDnum * chains,
	       SequencesReader *seqReadInfo, int WORDLENGTH)
{
	char *sequenceFilename = seqReadInfo->m_seqFilename;
	Annotation *annot = rdmaps->annotations;
	IDnum latestPreNodeID;
	InsertionMarker *currentMarker = insertionMarkers;
	IDnum sequenceIndex;
	Coordinate currentPosition, nextStop;
	IDnum preNodeCounter = 1;
	FILE *file = NULL;
	char line[50000];
	int lineLength = 50000;
	// Coordinate readIndex;
	boolean tooShort;
	Kmer initialKmer;
	char c;
	RoadMap *rdmap;
	IDnum annotIndex, lastAnnotIndex;
	IDnum markerIndex, lastMarkerIndex;

	char **sequences = malloc((sequenceCount_pg(preGraph) + 1) * sizeof(char *));
	sequenceIndex = 0;

	velvetLog("Sequence file = %s, is binary: %s\n", sequenceFilename, seqReadInfo->m_bIsBinary ? "true" : "false");
	if (!seqReadInfo->m_bIsBinary) {
		file = fopen(sequenceFilename, "r");
		if (file == NULL) 
			exitErrorf(EXIT_FAILURE, true, "Could not read %s", sequenceFilename);
		// Reading sequence descriptor in first line
		if (sequenceCount_pg(preGraph) > 0 && !fgets(line, lineLength, file))
			exitErrorf(EXIT_FAILURE, true, "%s incomplete.", sequenceFilename);
		seqReadInfo->m_pFile = file;

		char *currentSeq = NULL;
		size_t currentSeqSize = 0;
		while(fgets(line, lineLength, seqReadInfo->m_pFile) != NULL) {
			if (line[0] == '>') {
				// Start of a new sequence
				if (currentSeq != NULL) {
					// Store the completed sequence
					sequences[sequenceIndex++] = currentSeq;
					currentSeq = NULL;
					currentSeqSize = 0;
				}

				continue; // Skip the header line
			}

			// Remove newline characters
			line[strcspn(line, "\r\n")] = '\0';
			size_t lineLen = strlen(line);

			// Expand and append to current sequence buffer
			currentSeq = realloc(currentSeq, currentSeqSize + lineLen + 1);
			if (!currentSeq) {
				fprintf(stderr, "Error: Memory allocation failed\n");
				exit(1);
			}

			memcpy(currentSeq + currentSeqSize, line, lineLen);
			currentSeqSize += lineLen;
			currentSeq[currentSeqSize] = '\0';
		}

		// Add the last sequence
		if (currentSeq != NULL) {
			sequences[sequenceIndex++] = currentSeq;
		}

		// fclose(file);

		if (sequenceIndex != sequenceCount_pg(preGraph)) {
			fprintf(stderr, "Warning: Expected %d sequences, found %d\n", sequenceCount_pg(preGraph), sequenceIndex);
		}
	}
	velvetLog("Sequence 0 = %s\n", sequences[0]);
	velvetLog("Sequence 1 = %s\n", sequences[1]);
	velvetLog("Sequence %d = %s\n", sequenceCount_pg(preGraph) - 1, sequences[sequenceCount_pg(preGraph) - 1]);


	// Get starting annotation for all sequences with a prefix sum
	Annotation **startAnnot = malloc((sequenceCount_pg(preGraph) + 1) * sizeof(Annotation *));
	startAnnot[0] = annot;
	for (IDnum i = 1; i <= sequenceCount_pg(preGraph); i++) {
		startAnnot[i] = startAnnot[i - 1] + getAnnotationCount(getRoadMapInArray(rdmaps, i - 1));
	}

	// Get starting marker for all sequences with a prefix sum
	InsertionMarker **startMarkers = malloc((sequenceCount_pg(preGraph) + 1) * sizeof(InsertionMarker *));
	startMarkers[0] = currentMarker;
	for (IDnum i = 1; i <= sequenceCount_pg(preGraph); i++) {
		startMarkers[i] = startMarkers[i - 1] + markerCounters[i];
		// startMarkers[i] = startMarkers[i - 1] + markerCounters[i - 1];

	}

	#pragma omp parallel for schedule(dynamic) num_threads(2) private(rdmap, annotIndex, lastAnnotIndex, annot, markerIndex, lastMarkerIndex, currentPosition, currentMarker, latestPreNodeID, initialKmer, tooShort, nextStop)
	// Now that we have read all of the annotations, we go on to create the preNodes and tie them up
	for (sequenceIndex = 1; sequenceIndex <= sequenceCount_pg(preGraph); sequenceIndex++) {
		if (sequenceIndex % 1000000 == 1)
			velvetLog("Sequence %li / %li\n", (long) sequenceIndex,
			       (long) sequenceCount_pg(preGraph));
		if (sequenceIndex % 1000000 == 2)
		velvetLog("Sequence %li / %li\n", (long) sequenceIndex,
				(long) sequenceCount_pg(preGraph));

		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		annot = startAnnot[sequenceIndex - 1];
		annotIndex = 0;
		lastAnnotIndex = getAnnotationCount(rdmap);
		markerIndex = 0;
		lastMarkerIndex = markerCounters[sequenceIndex];
		currentPosition = 0;
		currentMarker = startMarkers[sequenceIndex - 1];
		IDnum threadLocalPreNodeCounter = chains[sequenceIndex];

		// Reading first (k-1) nucleotides
		tooShort = false;
		clearKmer(&initialKmer);
		//velvetLog("Initial kmer: ");
		TightString *tString = NULL;
		char *strString = NULL;
		if (seqReadInfo->m_bIsBinary) {
			tString = getTightStringInArray(seqReadInfo->m_sequences->tSequences, sequenceIndex - 1);
			strString = readTightString(tString);
		}

		Coordinate readIndex = 0;
		for (readIndex; readIndex < WORDLENGTH - 1; readIndex++) {
			if (seqReadInfo->m_bIsBinary) {
				if (readIndex >= tString->length) {
					tooShort = true;
					break;
				}

				c = strString[readIndex];
			} else {
				// Use sequences array for text input
				if (readIndex >= strlen(sequences[sequenceIndex - 1])) {
					tooShort = true;
					break;
				}
				c = sequences[sequenceIndex - 1][readIndex];
				switch (c) {
				case 'A':
				case 'N':
					pushNucleotide(&initialKmer, ADENINE);
					break;
				case 'C':
					pushNucleotide(&initialKmer, CYTOSINE);
					break;
				case 'G':
					pushNucleotide(&initialKmer, GUANINE);
					break;
				case 'T':
					pushNucleotide(&initialKmer, THYMINE);
					break;
				default:
					velvetLog
						("Irregular sequence file: are you sure your Sequence and Roadmap file come from the same source?\n");
					fflush(stdout);
					abort();
				}
			}
		}

		// Skipping sequence if too short
		if (tooShort) {
			//velvetLog("Skipping short read.. %d\n", sequenceIndex);
			// chains[sequenceIndex] = preNodeCounter;
			if (seqReadInfo->m_bIsBinary) {
				free(strString);
			} 
			// else {
			// 	if (!fgets(line, lineLength, file) && sequenceIndex < sequenceCount_pg(preGraph))
			// 		exitErrorf(EXIT_FAILURE, true, "%s incomplete.", sequenceFilename);
			// }
			continue;
		}

		// velvetLog("Sequence %li: past tooShort check, about to loop through annots\n", (long) sequenceIndex);
		char *currString = NULL;
		if (seqReadInfo->m_bIsBinary) {
			currString = &strString[readIndex];
			seqReadInfo->m_ppCurrString = &currString;
		}
		latestPreNodeID = 0;

		while (annotIndex < lastAnnotIndex) {
			if (markerIndex == lastMarkerIndex
			    || getPosition(annot) <=
			    getInsertionMarkerPosition(currentMarker))
				nextStop = getPosition(annot);
			else {
				nextStop = getInsertionMarkerPosition(currentMarker);
			}

			if (currentPosition != nextStop) {
				if (seqReadInfo->m_bIsBinary) {
					if (readIndex >= tString->length) {
						velvetLog("readIndex %ld beyond string len %ld\n", (uint64_t) readIndex, (uint64_t) tString->length);
						exit(1);
					}
				}
				
				// velvetLog("In annot loop, about to addPreNode, ID %d, readIndex %d, markerIndex %li, lastMarkerIndex %li, currentPosition %d, nextStop %d\n", threadLocalPreNodeCounter, readIndex, (long) markerIndex, (long) lastMarkerIndex, currentPosition, nextStop);
				addPreNodeToPreGraph_pgOMP(preGraph,
							currentPosition,
							nextStop,
							seqReadInfo,
							&initialKmer,
							threadLocalPreNodeCounter, sequences[sequenceIndex - 1], readIndex, sequenceIndex);
				threadLocalPreNodeCounter++;
				// if (latestPreNodeID == 0) {
				// 	chains[sequenceIndex] = threadLocalPreNodeCounter;
				// }
				latestPreNodeID = threadLocalPreNodeCounter;
				currentPosition = nextStop;
			}

			while (markerIndex < lastMarkerIndex
			       && getInsertionMarkerPosition(currentMarker)
			       == nextStop) {
				convertMarker(currentMarker,
					      latestPreNodeID);

				currentMarker++;
				markerIndex++;
			}

			// This part skips a bunch of nucleotides that are not important and moves the initialKmer window?
			Coordinate lastKmerIndex;
			while (annotIndex < lastAnnotIndex && getPosition(annot) == nextStop) {
				// Continue from lastKmerIndex
				lastKmerIndex = readIndex;
				for (readIndex; readIndex - lastKmerIndex < getAnnotationLength(annot); readIndex++) {
					if (seqReadInfo->m_bIsBinary) {
						c = *currString;
						currString += 1;   // increment the pointer
					} else {
						// Use sequences array for text input
						if (readIndex >= strlen(sequences[sequenceIndex - 1])) {
							velvetLog("readIndex %ld beyond string len %ld\n", (uint64_t) readIndex, (uint64_t) strlen(sequences[sequenceIndex - 1]));
							// exit(1);
							break;
						}
						c = sequences[sequenceIndex - 1][readIndex];

						switch (c) {
						case 'A':
						case 'N':
							pushNucleotide(&initialKmer, ADENINE);
							break;
						case 'C':
							pushNucleotide(&initialKmer, CYTOSINE);
							break;
						case 'G':
							pushNucleotide(&initialKmer, GUANINE);
							break;
						case 'T':
							pushNucleotide(&initialKmer, THYMINE);
							break;
						default:
							velvetLog
								("Irregular sequence file: are you sure your Sequence and Roadmap file come from the same source?\n");
							fflush(stdout);
#ifdef DEBUG 
							abort();
#endif 
							exit(1);
						}
					}
				}

				annot = getNextAnnotation(annot);
				annotIndex++;
			}

		}
		
		// velvetLog("Sequence %li: past annot loop, about to loop through rest of markers\n", (long) sequenceIndex);
		while (markerIndex < lastMarkerIndex) {
			// velvetLog("Sequence %li: In last marker while loop, markerIndex %li, lastMarkerIndex %li\n", (long) sequenceIndex, (long) markerIndex, (long) lastMarkerIndex);
			if (currentPosition ==
			    getInsertionMarkerPosition(currentMarker)) {
				convertMarker(currentMarker,
					      latestPreNodeID);
				currentMarker++;
				markerIndex++;
			} else {
				nextStop = getInsertionMarkerPosition(currentMarker);
				// velvetLog("Sequence %li: In last marker while loop, about to addPreNode, ID %d, readIndex %d, markerIndex %li, lastMarkerIndex %li, currentPosition %d, nextStop %d\n", (long) sequenceIndex, threadLocalPreNodeCounter, readIndex, (long) markerIndex, (long) lastMarkerIndex, currentPosition, nextStop);

				// threadLocalPreNodeCounter = __sync_fetch_and_add(&preNodeCounter, 1); 
				//	velvetLog("Adding pre nodes from %lli to %lli\n", (long long) currentPosition, (long long) nextStop);
				addPreNodeToPreGraph_pgOMP(preGraph,
							currentPosition,
							nextStop, seqReadInfo,
							&initialKmer,
							threadLocalPreNodeCounter, sequences[sequenceIndex - 1], readIndex, sequenceIndex);
				readIndex += nextStop - currentPosition;
				threadLocalPreNodeCounter++;
				// if (latestPreNodeID == 0)
				// 	chains[sequenceIndex] = threadLocalPreNodeCounter;
				latestPreNodeID = threadLocalPreNodeCounter;
				currentPosition = getInsertionMarkerPosition(currentMarker);
			}
		}

		// if (latestPreNodeID == 0)
		// 	chains[sequenceIndex] = threadLocalPreNodeCounter;
		
		// Check if threadLocalPredNodeCounter is equal to chains[sequenceIndex+1]
		if (sequenceIndex != sequenceCount_pg(preGraph) && threadLocalPreNodeCounter != chains[sequenceIndex + 1]) {
			velvetLog("ThreadLocalPreNodeCounter %li != chains[%li] %li\n", (long) threadLocalPreNodeCounter, (long) sequenceIndex + 1, (long) chains[sequenceIndex + 1]);
			exit(1);
		}
		// velvetLog("Sequence %li: Done\n", (long) sequenceIndex);
	}

	// preNodeCounter check
	velvetLog("preNodeCounter = %li, actual number of preNodes = %li\n", (long) preNodeCounter, (long) preNodeCount_pg(preGraph));
	// for (size_t i = 0; i < sequenceCount_pg(preGraph) + 1; i++) {
    //     free(sequences[i]);  // free each individual sequence string
    // }
	free(sequences);
	free(startAnnot);
	free(startMarkers); 
	free(markerCounters);
	if (!seqReadInfo->m_bIsBinary) {
		fclose(file);
	}

}

static void connectPreNodeToTheNext(IDnum * currentPreNodeID,
				    IDnum nextPreNodeID,
				    Coordinate * currentPosition,
				    IDnum sequenceIndex,
				    boolean isReference,
				    PreGraph * preGraph)
{
	if (nextPreNodeID == 0)
		return;

#ifdef _OPENMP
	lockTwoNodes(*currentPreNodeID, nextPreNodeID);
#endif

	if (isReference)
		incrementNodeReferenceMarkerCount_pg(preGraph, nextPreNodeID);

	if (!isReference && *currentPreNodeID != 0)
		createPreArc_pg(*currentPreNodeID, nextPreNodeID,
				preGraph);

#ifdef _OPENMP
	unLockTwoNodes(*currentPreNodeID, nextPreNodeID);
#endif

	*currentPreNodeID = nextPreNodeID;

	*currentPosition +=
	    getPreNodeLength_pg(*currentPreNodeID, preGraph);

}

static IDnum chooseNextInternalPreNode(IDnum currentPreNodeID,
				       IDnum sequenceIndex,
				       PreGraph * preGraph, IDnum * chains)
{
	if (currentPreNodeID >= preNodeCount_pg(preGraph))
		return 0;
	if (sequenceIndex >= sequenceCount_pg(preGraph))
		return currentPreNodeID + 1;
	if (currentPreNodeID + 1 < chains[sequenceIndex + 1])
		return currentPreNodeID + 1;
	return 0;
}

static void connectAnnotation(IDnum * currentPreNodeID, Annotation * annot,
			      Coordinate * currentPosition,
			      IDnum sequenceIndex, boolean isReference,
			      PreGraph * preGraph)
{
	IDnum nextPreNodeID = getStartID(annot);

	connectPreNodeToTheNext(currentPreNodeID, nextPreNodeID,
				currentPosition, 
				sequenceIndex, isReference, preGraph);

	while (*currentPreNodeID != getFinishID(annot)) {
		nextPreNodeID = (*currentPreNodeID) + 1;

		connectPreNodeToTheNext(currentPreNodeID, nextPreNodeID,
					currentPosition,
					sequenceIndex,
					isReference,
					preGraph);
	}
}

static void reConnectAnnotation(IDnum * currentPreNodeID, Annotation * annot,
			      Coordinate * currentPosition,
			      IDnum sequenceIndex, 
			      PreGraph * preGraph,
			      PreMarker ** previous)
{
	IDnum nextPreNodeID = getStartID(annot);

#ifdef _OPENMP
	lockNode(nextPreNodeID);
#endif
	*previous = addPreMarker_pg(preGraph, 
			nextPreNodeID,
			sequenceIndex,
			currentPosition, 
			*previous);
#ifdef _OPENMP
	unLockNode(nextPreNodeID);
#endif

	while (*currentPreNodeID != getFinishID(annot)) {
		nextPreNodeID = (*currentPreNodeID) + 1;

#ifdef _OPENMP
		lockNode(nextPreNodeID);
#endif
		*previous = addPreMarker_pg(preGraph, 
				nextPreNodeID,
				sequenceIndex,
				currentPosition,
				*previous);
#ifdef _OPENMP
		unLockNode(nextPreNodeID);
#endif
		*currentPreNodeID = nextPreNodeID;
	}
}

static void createPreMarkers(RoadMapArray * rdmaps, PreGraph * preGraph,
			    IDnum * chains)
{
	IDnum sequenceIndex;
	IDnum referenceCount = rdmaps->referenceCount;
#ifndef _OPENMP
	Annotation *annot = rdmaps->annotations;
#endif

#ifdef _OPENMP
	int threads = omp_get_max_threads();
	if (threads > 8)
		threads = 8;

	#pragma omp parallel for num_threads(threads)
#endif
	for (sequenceIndex = 1;
	     sequenceIndex <= referenceCount;
	     sequenceIndex++) {
#ifdef _OPENMP
		Annotation *annot = getAnnotationInArray(rdmaps->annotations, annotationOffset[sequenceIndex - 1]);
#endif
		RoadMap *rdmap;
		Coordinate currentPosition, currentInternalPosition;
		IDnum currentPreNodeID, nextInternalPreNodeID;
		IDnum annotIndex, lastAnnotIndex;
		PreMarker * previous;

		if (sequenceIndex % 1000000 == 0)
			velvetLog("Connecting %li / %li\n", (long) sequenceIndex,
			       (long) sequenceCount_pg(preGraph));

		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		annotIndex = 0;
		lastAnnotIndex = getAnnotationCount(rdmap);
		nextInternalPreNodeID = chooseNextInternalPreNode
		    (chains[sequenceIndex] - 1, sequenceIndex,
		     preGraph, chains);

		previous = NULL;
		currentPosition = 0;
		currentInternalPosition = 0;
		currentPreNodeID = 0;
		// Recursion up to last annotation
		while (annotIndex < lastAnnotIndex
		       || nextInternalPreNodeID != 0) {
			if (annotIndex == lastAnnotIndex
			    || (nextInternalPreNodeID != 0
				&& currentInternalPosition <
				getPosition(annot))) {
#ifdef _OPENMP
				lockNode(nextInternalPreNodeID);
#endif
				previous = addPreMarker_pg(preGraph, 
						nextInternalPreNodeID,
						sequenceIndex,
						&currentPosition,
						previous);
#ifdef _OPENMP
				unLockNode(nextInternalPreNodeID);
#endif
				currentPreNodeID = nextInternalPreNodeID;
				nextInternalPreNodeID =
				    chooseNextInternalPreNode
				    (currentPreNodeID, sequenceIndex,
				     preGraph, chains);
				currentInternalPosition +=
				    getPreNodeLength_pg(currentPreNodeID,
							preGraph);

			} else {
				reConnectAnnotation(&currentPreNodeID, annot,
						  &currentPosition,
						  sequenceIndex, 
						  preGraph,
						  &previous);
				annot = getNextAnnotation(annot);
				annotIndex++;
			}
		}
	}
}

// Threads each sequences and creates preArcs according to road map indications
static void connectPreNodes(RoadMapArray * rdmaps, PreGraph * preGraph,
			    IDnum * chains)
{
	IDnum sequenceIndex;
	IDnum referenceCount = rdmaps->referenceCount;
#ifdef _OPENMP
	annotationOffset = mallocOrExit(rdmaps->length + 1, Coordinate);
	annotationOffset[0] = 0;
	for (sequenceIndex = 1; sequenceIndex <= rdmaps->length; sequenceIndex++)
		annotationOffset[sequenceIndex] = annotationOffset[sequenceIndex - 1] +
						  getAnnotationCount(getRoadMapInArray(rdmaps, sequenceIndex - 1));
#else
	Annotation *annot = rdmaps->annotations;
#endif

	if (rdmaps->referenceCount > 0) 
		allocatePreMarkerCountSpace_pg(preGraph);

#ifdef _OPENMP
	int threads = omp_get_max_threads();
	if (threads > 8)
		threads = 8;

	#pragma omp parallel for num_threads(threads)
#endif
	for (sequenceIndex = 1;
	     sequenceIndex <= sequenceCount_pg(preGraph);
	     sequenceIndex++) {
#ifdef _OPENMP
		Annotation *annot = getAnnotationInArray(rdmaps->annotations, annotationOffset[sequenceIndex - 1]);
#endif
		RoadMap *rdmap;
		Coordinate currentPosition, currentInternalPosition;
		IDnum currentPreNodeID, nextInternalPreNodeID;
		IDnum annotIndex, lastAnnotIndex;
		boolean isReference;

		if (sequenceIndex % 1000000 == 0)
			velvetLog("Connecting %li / %li\n", (long) sequenceIndex,
			       (long) sequenceCount_pg(preGraph));

		rdmap = getRoadMapInArray(rdmaps, sequenceIndex - 1);
		annotIndex = 0;
		lastAnnotIndex = getAnnotationCount(rdmap);
		nextInternalPreNodeID = chooseNextInternalPreNode
		    (chains[sequenceIndex] - 1, sequenceIndex,
		     preGraph, chains);
		isReference = (sequenceIndex <= referenceCount);

		currentPosition = 0;
		currentInternalPosition = 0;
		currentPreNodeID = 0;
		// Recursion up to last annotation
		while (annotIndex < lastAnnotIndex
		       || nextInternalPreNodeID != 0) {
			if (annotIndex == lastAnnotIndex
			    || (nextInternalPreNodeID != 0
				&& currentInternalPosition <
				getPosition(annot))) {
				connectPreNodeToTheNext(&currentPreNodeID,
							nextInternalPreNodeID,
							&currentPosition,
							sequenceIndex,
							isReference,
							preGraph);
				nextInternalPreNodeID =
				    chooseNextInternalPreNode
				    (currentPreNodeID, sequenceIndex,
				     preGraph, chains);
				currentInternalPosition +=
				    getPreNodeLength_pg(currentPreNodeID,
							preGraph);

			} else {
				connectAnnotation(&currentPreNodeID, annot,
						  &currentPosition,
						  sequenceIndex, isReference,
						  preGraph);
				annot = getNextAnnotation(annot);
				annotIndex++;
			}
		}
	}

	if (rdmaps->referenceCount > 0) {
		allocatePreMarkerSpace_pg(preGraph);
		createPreMarkers(rdmaps, preGraph, chains);	
	}

#ifdef _OPENMP
	free(annotationOffset);
	annotationOffset = NULL;
#endif
}

// Post construction memory deallocation routine (of sorts, could certainly be optimized)
static void
cleanUpMemory(PreGraph * preGraph, RoadMapArray * rdmaps, IDnum * chains)
{
	// Killing off roadmaps
	destroyRoadMapArray(rdmaps);

	// Finishing off the chain markers
	free(chains);
}

// The full monty, wrapped up in one function
PreGraph *newPreGraph_pg(RoadMapArray * rdmapArray, SequencesReader *seqReadInfo)
{
	int WORDLENGTH = rdmapArray->WORDLENGTH;
	IDnum sequenceCount = rdmapArray->length;
	IDnum *markerCounters = callocOrExit(sequenceCount + 1, IDnum);
	IDnum *chains = callocOrExit(sequenceCount + 1, IDnum);
	InsertionMarker *insertionMarkers;
	InsertionMarker *veryLastMarker;

	PreGraph *preGraph =
	    emptyPreGraph_pg(sequenceCount, rdmapArray->referenceCount, rdmapArray->WORDLENGTH, rdmapArray->double_strand);

	velvetLog("Creating insertion markers\n");
	setInsertionMarkers(rdmapArray, markerCounters, &veryLastMarker,
			    &insertionMarkers);

	velvetLog("Counting preNodes\n");
#ifdef _OPENMP
	countPreNodesOMP(rdmapArray, preGraph, markerCounters,
		      insertionMarkers, veryLastMarker, chains);
#else
	countPreNodes(rdmapArray, preGraph, markerCounters,
		      insertionMarkers, veryLastMarker);
#endif

	velvetLog("%li preNodes counted, creating them now\n",
	       (long) preNodeCount_pg(preGraph));
#ifdef _OPENMP
	createPreNodesOMP(rdmapArray, preGraph, markerCounters,
		       insertionMarkers, veryLastMarker, chains,
		       seqReadInfo, WORDLENGTH);
#else
	createPreNodes(rdmapArray, preGraph, markerCounters,
		       insertionMarkers, veryLastMarker, chains,
		       seqReadInfo, WORDLENGTH);
#endif

	velvetLog("Adjusting marker info...\n");
	convertInsertionMarkers(insertionMarkers, veryLastMarker, chains);

#ifdef _OPENMP
	createNodeLocks(preGraph);
#endif
	velvetLog("Connecting preNodes\n");
	connectPreNodes(rdmapArray, preGraph, chains);

	velvetLog("Cleaning up memory\n");
	cleanUpMemory(preGraph, rdmapArray, chains);
#ifdef _OPENMP
	free(nodeLocks);
	nodeLocks = NULL;
#endif 

	velvetLog("Done creating preGraph\n");

	return preGraph;
}
