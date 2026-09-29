# Why run locally?

## Cost Comparison

| Approach | Cost | What You Get | Data Privacy |
|---|---|---|---|
| **This pipeline** | $0 (free, open source) | 33 default + 1 opt-in analysis steps | Your data never leaves your machine |
| Clinical WGS interpretation | $500-5,000 | 1-page report, selected genes only | Lab retains your data |
| Nebula/Dante report | $0-200 (included/add-on) | Web dashboard, limited depth | Data on company servers |
| 23andMe Health | $229 | ~10 health reports from array data | Data shared with research partners |
| Genetic counselor consultation | $200-500/hour | Expert interpretation of specific findings | HIPAA-protected |

**The pipeline is complementary, not a replacement.** Use it for comprehensive self-analysis, then bring specific findings to a genetic counselor or physician for clinical interpretation.

## Privacy and Security

Your genome is the most permanent piece of personal data you have. Unlike a password, you cannot change it if it leaks.

**This pipeline keeps your data local:**
- No data is uploaded to any server, cloud, or API
- No telemetry, no analytics, no tracking
- After completing [reference setup](00-reference-setup.md) (which pre-downloads all databases and Docker images), the core pipeline runs offline
- A few steps fetch small public resources on first use if not already present: step 4b downloads the ENCODE blacklist (~50 KB), step 13 can download the VEP cache (~26 GB) if skipped during setup, and steps 25/26 download scoring files and reference panels from public FTP servers. No sample data is ever uploaded. All downloads are cached after the first run

**Recommendations for securing your data:**
- Store genomic data on an encrypted filesystem (LUKS on Linux, FileVault on macOS, BitLocker on Windows)
- Never upload raw FASTQ/BAM/VCF files to unencrypted cloud storage
- If using a NAS, enable encryption at rest
- Be cautious with VCF files in particular — they are small enough to accidentally email or upload
- Consider who has physical access to the machine where your data is stored
- If you delete your data, use `shred` (Linux) or secure erase — standard file deletion leaves data recoverable

**GDPR note:** If you are in the EU, your genomic data is classified as "special category personal data" under GDPR Article 9. Processing it locally for personal use is lawful. Sharing it with third parties (including cloud services) may require explicit consent and appropriate safeguards.

