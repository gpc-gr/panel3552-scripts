package org.analysis;

import java.io.File;
import java.io.PrintStream;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class DuplicatedReadRemover {

    public static void main(String[] args) {
        String bamFile = args[0];
        new DuplicatedReadRemover().removeDuplicatedReads(bamFile, System.out);
    }

    public void removeDuplicatedReads(String bamFile, PrintStream printStream) {

        SAMFileReader samFileReader= new SAMFileReader(new File(bamFile));
        samFileReader.setValidationStringency(ValidationStringency.SILENT);

        SAMFileWriter samFileWriter = new SAMFileWriter(samFileReader.getFileHeader(), printStream);
        SAMRecordIterator samRecordIterator = samFileReader.iterator();

        String currentReadName = null;
        int count = 0;

        double[] rates       = new double[]{0, 0};
        SAMRecord[] readPair = new SAMRecord[]{null, null};

        while (samRecordIterator.hasNext()) {
            count++;
            SAMRecord samRecord = samRecordIterator.next();
            if ((count % 1000000) == 0) {
                System.err.println(count + " records have been processed");
            }
            if (isEmptyPair(readPair)) {
                currentReadName = samRecord.getReadName();
                addRead(samRecord, readPair, rates);
            }
            else if (currentReadName.equals(samRecord.getReadName())) {
                addRead(samRecord, readPair, rates);
            }
            else {
                writeSequence(readPair, samFileWriter);
                rates[0]    = rates[1]    = 0;
                readPair[0] = readPair[1] = null;
                currentReadName = samRecord.getReadName();
                addRead(samRecord, readPair, rates);
            }
        }
        if (isEmptyPair(readPair) == false) {
            writeSequence(readPair, samFileWriter);
        }
        samFileReader.close();
    }

    private boolean isEmptyPair(SAMRecord[] readPair) {
        if (readPair[0] != null) {
            return false;
        }
        if (readPair[1] != null) {
            return false;
        }
        return true;
    }

    private void addRead(SAMRecord samRecord, SAMRecord[] readPair, double[] rates) {
        double rate = getReadRate(samRecord);
        boolean firstOfPairFlag = samRecord.getFirstOfPairFlag();
        if (firstOfPairFlag) {
            if (readPair[0] == null || rate > rates[0]) {
                rates[0]    = rate;
                readPair[0] = samRecord;
            }
        }
        else {
            if (readPair[1] == null || rate > rates[1]) {
                rates[1]    = rate;
                readPair[1] = samRecord;
            }
        }
    }

    private void writeSequence(SAMRecord[] readPair, SAMFileWriter samFileWriter) {
        if (readPair[0] != null) {
            samFileWriter.addAlignment(readPair[0]);
        }
        if (readPair[1] != null) {
            samFileWriter.addAlignment(readPair[1]);
        }
    }

    private double getReadRate(SAMRecord samRecord) {
        if (samRecord.getReadString().equals("*")) {
            return 0;
        }
        Cigar cigar = samRecord.getCigar();
        int requiredBaseCount = 0;

        for (CigarElement cigarElement : cigar.getCigarElements()) {
            if (cigarElement.getOperator() == CigarOperator.MATCH_OR_MISMATCH) {
                requiredBaseCount += cigarElement.getLength();
            }
            else if (cigarElement.getOperator() == CigarOperator.INSERTION) {
                requiredBaseCount += cigarElement.getLength();
            }
            else if (cigarElement.getOperator() == CigarOperator.SOFT_CLIP) {
                requiredBaseCount += cigarElement.getLength();
            }
            else if (cigarElement.getOperator() == CigarOperator.HARD_CLIP) {
                requiredBaseCount += cigarElement.getLength();
            }
        }
        return samRecord.getReadBases().length / ((double)requiredBaseCount);
    }
}
