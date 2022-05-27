.. automodule:: scReadSim

API
===

Import scReadSim as::

   import scReadSim

Utility
~~~~~~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: _autosummary

   scReadSim.scReadSim.Utility.CallPeak
   scReadSim.scReadSim.Utility.ExtractBAMCoverage
   scReadSim.scReadSim.Utility.scATAC_CreateFeatureSets
   scReadSim.scReadSim.Utility.scRNA_CreateFeatureSets
   scReadSim.scReadSim.Utility.bam2countmat
   scReadSim.scReadSim.Utility.bam2countmat_INPUT
   scReadSim.scReadSim.Utility.find_nearest
   scReadSim.scReadSim.Utility.match_peak
   scReadSim.scReadSim.Utility.ComplementFeature


scATAC_GenerateBAM
~~~~~~~~~~~~~~~~~~
.. autosummary::
   :toctree: _autosummary

   scReadSim.scReadSim.scATAC_GenerateBAM.flatten
   scReadSim.scReadSim.scATAC_GenerateBAM.cellbarode_generator
   scReadSim.scReadSim.scATAC_GenerateBAM.scATAC_INPUT_PerTruePeakEdition
   scReadSim.scReadSim.scATAC_GenerateBAM.scATAC_PerTruePeakEdition
   scReadSim.scReadSim.scATAC_GenerateBAM.scATAC_SampleSyntheticReads
   scReadSim.scReadSim.scATAC_GenerateBAM.scATAC_GenerateBAMCoord
   scReadSim.scReadSim.scATAC_GenerateBAM.scATAC_SampleSyntheticReads_INPUT
   scReadSim.scReadSim.scATAC_GenerateBAM.scATAC_GenerateBAMCoord_INPUT
   scReadSim.scReadSim.scATAC_GenerateBAM.scATAC_CombineBED
   scReadSim.scReadSim.scATAC_GenerateBAM.scATAC_BED2FASTQ
   scReadSim.scReadSim.scATAC_GenerateBAM.AlignSyntheticBam_Pair


scRNA_GenerateBAM
~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary

   scReadSim.scReadSim.scRNA_GenerateBAM.flatten
   scReadSim.scReadSim.scRNA_GenerateBAM.cellbarode_generator
   scReadSim.scReadSim.scRNA_GenerateBAM.scRNA_SampleSyntheticReads
   scReadSim.scReadSim.scRNA_GenerateBAM.scRNA_PerTruePeakEdition
   scReadSim.scReadSim.scRNA_GenerateBAM.scRNA_GenerateBAMCoord
   scReadSim.scReadSim.scRNA_GenerateBAM.scRNA_CombineBED
   scReadSim.scReadSim.scRNA_GenerateBAM.scRNA_BED2FASTQ
   scReadSim.scReadSim.scRNA_GenerateBAM.AlignSyntheticBam_Pair

GenerateSyntheticCount
~~~~~~~~~~~~~~~~~~~~~~

.. autosummary::
   :toctree: _autosummary

   scReadSim.scReadSim.GenerateSyntheticCount.scATAC_GenerateSyntheticCount
   scReadSim.scReadSim.GenerateSyntheticCount.scRNA_GenerateSyntheticCount
