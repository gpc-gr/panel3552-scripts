package org.analysis;

import java.io.File;
import java.io.PrintStream;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;



public class SAMFileWriter {

    private PrintStream printStream;

    public SAMFileWriter(SAMFileHeader samFileHeader, PrintStream printStream) {
        this.printStream = printStream;
        setSAMFileHeader(samFileHeader);
    }

    public SAMFileWriter(PrintStream printStream) {
        this.printStream = printStream;
    }

    public void setSAMFileHeader(SAMFileHeader samFileHeader) {
        this.printStream.print(samFileHeader.getTextHeader());
    }

    public void addAlignment(SAMRecord samRecord) {
        this.printStream.print(samRecord.getSAMString());
    }
}
