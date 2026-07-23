# Reference Annotation Cache guide

The cache stores assembly accessions in `GenomeCache.sqlite` and downloaded GFF annotations in `Cache/GFF/`. It does not store genome FASTA sequences. Use this application to inspect references, verify integrity, export a portable ZIP, import another cache and merge compatible records. Create a backup before destructive cleanup and prefer portable ZIP packages over isolated SQLite files.
