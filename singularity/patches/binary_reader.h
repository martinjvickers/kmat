/*
This file is a part of KMC software distributed under GNU GPL 3 licence.
The homepage of the KMC project is http://sun.aei.polsl.pl/kmc

Authors: Sebastian Deorowicz, Agnieszka Debudaj-Grabysz, Marek Kokot

Version: 3.2.4
Date   : 2024-02-09
*/

#ifndef _BINARY_READER_H
#define _BINARY_READER_H

#include "defs.h"
#include "params.h"
#include "queues.h"
#include "percent_progress.h"
#define DONT_DEFINE_UCHAR
#include "../kmc_api/kmc_file.h"
#include "critical_error_handler.h"
#include "bam_utils.h"
#include <sys/stat.h>
#include <zlib.h>

// kmat patch: decompress .gz with zlib gzFile (same path as zcat) and feed plain
// FASTQ packs to readers. Avoids KMC's manual inflate loop + Cloudflare zlib,
// which fail on many real Illumina / pigz multi-member .fq.gz files with:
//   "Some error while reading gzip file ... fastq_reader.cpp"

class CBinaryFilesReader
{
	bool is_file(const char* path)
	{
#ifdef _WIN32
		typedef struct _stat64 stat_struct;
		const auto& stat_func = _stat64;
#else
		typedef struct stat stat_struct;
		const auto& stat_func = stat;
#endif
		stat_struct buf;
		if (stat_func(path, &buf) == -1)
			return false;
		return (buf.st_mode & S_IFMT) == S_IFREG;
	}
	uint32 part_size;
	CInputFilesQueue* input_files_queue;
	CMemoryPool *pmm_binary_file_reader;
	vector<CBinaryPackQueue*> binary_pack_queues;
	CBamTaskManager* bam_task_manager = nullptr;
	uint64 total_size;
	uint64 predicted_size;	
	CPercentProgress percent_progress;

	InputType input_type; //for bam input behaviour of this class is quite different, for example only one file is readed at once
						  //also for KMC, where this class does almost nothing, just sends kmc file path to reader
	void notify_readed(uint64 readed)
	{		
		percent_progress.NotifyProgress(readed);
	}
	CompressionType get_compression_type(const string& name)
	{
		if (name.size() > 3 && string(name.end() - 3, name.end()) == ".gz")
			return CompressionType::gzip;
		else
			return CompressionType::plain;
	}

	struct OpenInput
	{
		FILE* fp = nullptr;
		gzFile gzf = nullptr;
		CBinaryPackQueue* q = nullptr;
		// Mode pushed into the pack queue. .gz is pre-decompressed via gzread,
		// so readers always see CompressionType::plain.
		CompressionType pack_mode = CompressionType::plain;
	};

	void CloseInput(OpenInput& in)
	{
		if (in.gzf)
		{
			gzclose(in.gzf);
			in.gzf = nullptr;
		}
		if (in.fp)
		{
			fclose(in.fp);
			in.fp = nullptr;
		}
	}

	void OpenFile(const string& file_name, OpenInput& in)
	{
		CloseInput(in);
		const CompressionType file_comp = get_compression_type(file_name);
		if (file_comp == CompressionType::gzip)
		{
			in.gzf = gzopen(file_name.c_str(), "rb");
			if (!in.gzf)
			{
				std::ostringstream ostr;
				ostr << "Error: cannot open gzip file: " << file_name << " for reading";
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			gzbuffer(in.gzf, 1 << 20);
			in.pack_mode = CompressionType::plain;
			return;
		}

		in.fp = fopen(file_name.c_str(), "rb");
		if (!in.fp)
		{
			std::ostringstream ostr;
			ostr << "Error: cannot open file: " << file_name << " for reading";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		setvbuf(in.fp, nullptr, _IONBF, 0);
		in.pack_mode = CompressionType::plain;
	}

	uint64 ReadInput(OpenInput& in, uchar* part, uint32 part_size)
	{
		if (in.gzf)
		{
			const int n = gzread(in.gzf, part, part_size);
			if (n < 0)
			{
				int err = Z_OK;
				const char* msg = gzerror(in.gzf, &err);
				std::ostringstream ostr;
				ostr << "Error while reading gzip file via zlib gzread";
				if (msg)
					ostr << ": " << msg;
				CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
			}
			return static_cast<uint64>(n);
		}
		return fread(part, 1, part_size, in.fp);
	}

	uint64_t skipSingleBGZFBlock(uchar* buff)
	{
		uint64_t pos = 0;
		if (buff[0] == 0x1f && buff[1] == 0x8b)
			;
		else
		{
			std::ostringstream ostr;
			ostr << "Fail: this is not gzip file";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}

		if (buff[2] != 8)
		{
			std::ostringstream ostr;
			ostr << "Error: CM flag is set to " << buff[2] << " instead of 8";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}

		if (!((buff[3] >> 2) & 1))
		{
			std::ostringstream ostr;
			ostr << "Error: FLG.FEXTRA is not set";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}

		pos = 10;

		uint16_t XLEN;
		read_uint16_t(XLEN, buff, pos);

		uchar SI1 = buff[pos++];
		uchar SI2 = buff[pos++];
		if (SI1 != 66 || SI2 != 67)
		{
			std::ostringstream ostr;
			ostr << "Error: SI1 != 66 || SI2 != 67";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		uint16_t LEN;
		read_uint16_t(LEN, buff, pos);
		if (LEN != 2)
		{
			std::ostringstream ostr;
			ostr << "Error: SLEN is " << LEN << " instead of 2";
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}

		uint16_t BGZF_block_size_tmp;
		read_uint16_t(BGZF_block_size_tmp, buff, pos);

		uint64_t BGZF_block_size = BGZF_block_size_tmp + 1ull;  //This integer gives the size of the containing BGZF block minus one.
		return BGZF_block_size;
	}

	uint64_t findLastBGZFBlockEnd(uchar* data, uint64_t size)
	{
		//header of BGZF block is 18 Bytes long
		uint64_t pos = 0;
		while (pos + 18 < size)
		{
			uint64_t npos = skipSingleBGZFBlock(data + pos);
			if (pos + npos > size)
				break;
			pos += npos;
		}
		return pos;
	}

	void ProcessSingleBamFile(const string& fname, uint32 file_no, uint32& id, bool& forced_to_finish)
	{
		FILE* file = fopen(fname.c_str(), "rb");
		if (!file)
		{
			std::ostringstream ostr;
			ostr << "Error: cannot open file " << fname;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		setvbuf(file, nullptr, _IONBF, 0);

		unsigned char eof_marker[] = { 0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

		unsigned char eof_to_chech[sizeof(eof_marker)];
		fseek(file, -static_cast<int>(sizeof(eof_marker)), SEEK_END);
		if (sizeof(eof_marker) != fread(eof_to_chech, 1, sizeof(eof_marker), file))
		{
			std::ostringstream ostr;
			ostr << "Error: cannot check EOF marker of BAM file: " << fname;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		if (!equal(begin(eof_marker), end(eof_marker), begin(eof_to_chech)))
		{
			std::ostringstream ostr;
			ostr << "Error: wrong EOF marker of BAM file: " << fname;
			CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
		}
		fseek(file, 0, SEEK_SET);

		uchar* data;
		pmm_binary_file_reader->reserve(data);

		uint64 size = 0;
		
		uint64 readed;
		
		while (!forced_to_finish && (readed = fread(data + size, 1, part_size - size, file)))
		{
			notify_readed(readed);
			size += readed;
			uint64_t lastBGFBlockEnd = findLastBGZFBlockEnd(data, size);
			uchar* newData;
			pmm_binary_file_reader->reserve(newData);

			uint64_t tail = size - lastBGFBlockEnd;
			memcpy(newData, data + lastBGFBlockEnd, tail);
			size = lastBGFBlockEnd;
			
			if (!bam_task_manager->PushBinaryPack(data, size, id, file_no))
			{
				pmm_binary_file_reader->free(data);
				forced_to_finish = true;
			}
			else
				id++;
			
			data = newData;
			size = tail;
		}
		if (!bam_task_manager->PushBinaryPack(data, size, id, file_no)) //last, possibly empty
		{
			pmm_binary_file_reader->free(data);			
			forced_to_finish = true;
		}
		else
			id++;
		fclose(file);
	}

	void ProcessBam()
	{
		uint32 file_no = 0;
		uint32 id = 0;
		string fname;
		bool forced_to_finish = false;
		
		while (!forced_to_finish && input_files_queue->pop(fname))
		{
			ProcessSingleBamFile(fname, file_no, id, forced_to_finish);
			++file_no;
		}
		bam_task_manager->NotifyBinaryReaderCompleted(id-1);		
	}

	/*
	* TODO: for now very simple version, that just reads k-mers from databases and transforms it to pseudoreads
	* it uses KMC API. Better solution would be to read binary data from kmc pre and suf files and pass it to the reader that would decode 
	* and prepare for splitter.
	* TODO: also this thread is now much more busy and probably should be counted to the numbre of threads
	*/
	void ProcessKMC()
	{
		std::string file_name;
		//for now im using only one binary_pack_queue, but it is just first and simple implementation
		notify_readed(0);
		while (input_files_queue->pop(file_name))
		{
			uchar* part;
			pmm_binary_file_reader->reserve(part);
			CKMCFile kmc_file;
			kmc_file.OpenForListing(file_name);
			CKMCFileInfo info;
			kmc_file.Info(info);
			uint32_t rec_size = 3 + info.kmer_length; // `>` + '\n' + kmer + '\n'
			uint32 kmers_in_part = part_size / rec_size;

			CKmerAPI kmer(info.kmer_length);

			uint32 id{};
			FilePart file_part = FilePart::Begin;
			uint32 c;
			uchar* ptr = part;
			bool forced_to_finish = false;
			while (kmc_file.ReadNextKmer(kmer, c))
			{
				*ptr++ = '>';
				*ptr++ = '\n';
				kmer.to_string((char*)ptr);
				ptr += info.kmer_length;
				*ptr++ = '\n';
				++id;
				if (id == kmers_in_part)
				{
					notify_readed(id);
					if (!binary_pack_queues[0]->push(part, id * rec_size, file_part, CompressionType::plain))
					{
						forced_to_finish = true;
						pmm_binary_file_reader->free(part);
						break;
					}
					id = 0;
					file_part = FilePart::Middle;
					pmm_binary_file_reader->reserve(part);
					ptr = part;
				}
			}

			if (id && !forced_to_finish)
			{
				notify_readed(id);
				if (!binary_pack_queues[0]->push(part, id * rec_size, file_part, CompressionType::plain))
				{
					pmm_binary_file_reader->free(part);
				}
				//id = 0;
				//file_part = FilePart::Middle;
				//pmm_binary_file_reader->reserve(part);
				//ptr = part;
			}
		}

		binary_pack_queues[0]->push(nullptr, 0, FilePart::End, CompressionType::plain);

		//user may specify more fastq_readers than input files
		for (auto& e : binary_pack_queues)
			e->mark_completed();
	}


public:
	CBinaryFilesReader(CKMCParams &Params, CKMCQueues &Queues, bool _show_progress)
		:
		percent_progress("Stage 1: ", _show_progress, Params.percentProgressObserver)
	{
		part_size = (uint32)Params.mem_part_pmm_binary_file_reader;
		input_files_queue = Queues.input_files_queue.get();
		pmm_binary_file_reader = Queues.pmm_binary_file_reader.get();
		//binary_pack_queues = Queues.binary_pack_queues;
		std::transform(
			Queues.binary_pack_queues.begin(),
			Queues.binary_pack_queues.end(),
			std::back_inserter(binary_pack_queues),
			[](const auto& uptr) {return uptr.get(); });

		bam_task_manager = Queues.bam_task_manager.get();
		auto files_copy = input_files_queue->GetCopy();		
		total_size = 0;
		predicted_size = 0;
		input_type = Params.file_type;

		while (!files_copy.empty())
		{
			string& f_name = files_copy.front();
			uint64 fsize;
			if(input_type == InputType::KMC)
			{
				CKMCFile kmc_file;
				if (!kmc_file.OpenForListing(f_name))
				{
					std::ostringstream ostr;
					ostr << "Cannot open KMC database: " << f_name;
					CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
				}				
				CKMCFileInfo info;
				kmc_file.Info(info);
				fsize = info.total_kmers;
				kmc_file.Close();
			}
			else
			{
				FILE* f = fopen(f_name.c_str(), "rb");
				if (!is_file(f_name.c_str()))
				{
					std::ostringstream ostr;
					ostr << "Error: " << f_name << " is not a file";
					CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
				}
				if (!f)
				{
					std::ostringstream ostr;
					ostr << "Cannot open file: " << f_name;
					CCriticalErrorHandler::Inst().HandleCriticalError(ostr.str());
				}
				my_fseek(f, 0, SEEK_END);
				fsize = my_ftell(f);
				fclose(f);
			}

			total_size += fsize;

			if (input_type == InputType::BAM)
			{
				predicted_size += (uint64)(0.7 * fsize); //TODO: is this correct?
			}
			else if (input_type == InputType::KMC)
			{
				predicted_size = fsize;
			}
			else
			{
				CompressionType compression = get_compression_type(f_name);
				switch (compression)
				{
				case CompressionType::plain:
					predicted_size += fsize;
					break;
				case CompressionType::gzip:
					predicted_size += (uint64)(3.2 * fsize);
					break;
				default:
					break;
				}
			}
			files_copy.pop();
		}
		percent_progress.SetMaxVal(total_size);
	}
	uint64 GetTotalSize()
	{
		return total_size;
	}
	uint64 GetPredictedSize()
	{
		return predicted_size;
	}

	void Process()
	{
		if (input_type == InputType::BAM)
		{
			ProcessBam();
			return;
		}
		if (input_type == InputType::KMC)
		{
			ProcessKMC();
			return;
		}
		std::string file_name;
		vector<OpenInput> files;
		files.reserve(binary_pack_queues.size());
		uchar* part = nullptr;

		uint32 completed = 0;
		notify_readed(0);

		for (uint32 i = 0; i < binary_pack_queues.size() && input_files_queue->pop(file_name); ++i)
		{
			OpenInput in;
			in.q = binary_pack_queues[i];
			OpenFile(file_name, in);
			files.push_back(in);
			pmm_binary_file_reader->reserve(part);
			uint64 readed = ReadInput(files.back(), part, part_size);
			notify_readed(readed);
			if (!files.back().q->push(part, readed, FilePart::Begin, files.back().pack_mode))
			{
				pmm_binary_file_reader->free(part);
				CloseInput(files.back());
				++completed;
			}
		}

		bool forced_to_finish = false;

		while (completed < files.size() && !forced_to_finish)
		{
			for (auto& f : files)
			{
				if (!f.fp && !f.gzf)
					continue;

				pmm_binary_file_reader->reserve(part);
				uint64 readed = ReadInput(f, part, part_size);
				notify_readed(readed);
				if (readed == 0) //end of file, need to open next one if exists
				{
					pmm_binary_file_reader->free(part);
					if (!f.q->push(nullptr, 0, FilePart::End, f.pack_mode))
					{
						forced_to_finish = true;
						break;
					}
					CloseInput(f);

					if (input_files_queue->pop(file_name))
					{
						OpenFile(file_name, f);
						pmm_binary_file_reader->reserve(part);
						readed = ReadInput(f, part, part_size);
						notify_readed(readed);
						if (!f.q->push(part, readed, FilePart::Begin, f.pack_mode))
						{
							pmm_binary_file_reader->free(part);
							forced_to_finish = true;
							break;
						}
					}
					else
					{
						++completed;
						f.q->mark_completed();
					}
				}
				else
				{
					if (!f.q->push(part, readed, FilePart::Middle, f.pack_mode))
					{
						pmm_binary_file_reader->free(part);
						forced_to_finish = true;
						break;
					}
				}
			}
		}

		for (auto& f : files)
			CloseInput(f);

		//user may specify more fastq_readers than input files
		for (auto& e : binary_pack_queues)
			e->mark_completed();
	}
};


class CWBinaryFilesReader
{
	std::unique_ptr<CBinaryFilesReader> reader;
public:
	CWBinaryFilesReader(CKMCParams &Params, CKMCQueues &Queues, bool show_progress = true)
	{
		reader = std::make_unique<CBinaryFilesReader>(Params, Queues, show_progress);
	}

	uint64 GetPredictedSize()
	{
		return reader->GetPredictedSize();
	}
	uint64 GetTotalSize()
	{
		return reader->GetTotalSize();
	}

	void operator()()
	{
		reader->Process();
	}
};

#endif

// ***** EOF