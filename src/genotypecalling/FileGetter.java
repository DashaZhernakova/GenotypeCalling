package genotypecalling;


import umcg.genetica.io.text.TextFile;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

public class FileGetter {
    public String[] fileNames;
    public String[] sampleNames;

    /**
     * Gets paths to files if a list of samples and corresponding file paths is given (file paths from the 2nd column)
     * @param fileListPath - path to the file with sample names in the first column and paths to snvmix files in the second column (tab as a separator)
     * @throws java.io.IOException
     */
    public void loadFileNames(String fileListPath) throws IOException {
        System.out.println("Getting files from file list: " + fileListPath);
        TextFile filelist = new TextFile(fileListPath, false);
        int col = 1;
        if (filelist.countCols(TextFile.tab) == 1){
            col = 0;
        }
        System.out.println("Getting file names from column " + col);

        fileNames = filelist.readAsArray(col, TextFile.tab);
        System.out.println("Found " + fileNames.length + " files");
        filelist.close();
    }

    /**
     * Gets files from folder dirName whose names match snvmixFnamePattern (regex)
     *
     * @param dirName
     * @param genPattern regexp pattern for matching genotype files
     *
     */
    public void loadFileNames(String dirName, String genPattern){
        System.out.println("Getting files from: " + dirName + " \tmatching pattern: " + genPattern);

        File dir = new File(dirName);
        ArrayList<String> fnames = new ArrayList<String>();
        for (File ch : dir.listFiles())
            if (ch.isDirectory())
                for (File child : ch.listFiles())
                    if ((child.getName().matches(genPattern)) && (child.length() != 0))
                        fnames.add(child.getPath());
        fileNames = fnames.toArray(new String[0]);
        System.out.println("Found " + fileNames.length + " files");

    }

    /**
     * Gets genotype files from dirName, matching genPattern, ordered by chromosome
     * @param dirName
     * @param genPattern
     */
    public void loadFileNamesByChr(String dirName, String genPattern){
        System.out.println("Getting files from: " + dirName + " \tmatching pattern: " + genPattern);
        System.out.println("Files will be sorted by chromosome");
        File dir = new File(dirName);
        String chrString;
        String pattern = null;
        fileNames = new String[24];
        for (int chr = 1; chr < 24; chr++){
            if (chr == 23)
                chrString = "X";
            else
                chrString = String.valueOf(chr);
            System.out.println("Chromosome " + chrString);
            pattern = genPattern.replace("[1-9X]+", chrString);

            if (pattern == null)
                pattern = genPattern.replace("(chr)", chrString).replace(".", "\\.");
            if (pattern == null){
                System.out.println("Wrong pattern: " + genPattern);
                System.exit(-1);
            }
            System.out.println("Pattern: " + pattern);
            for (File f : dir.listFiles()) {
                if (f.getName().matches(pattern)) {
                    fileNames[chr] = f.getPath();

                }
            }
        }
        System.out.println("Found " + fileNames.length + " files");

    }

    /**
     * Gets sample names if a list of samples and corresponding file paths is given (sample names from the 1st column)
     * @param fileListPath - path to the file with sample names in the first column and paths to snvmix files in the second column (tab as a separator)
     * @throws java.io.IOException
     */
    public void loadSampleNames(String fileListPath) throws IOException {
        System.out.println("Getting sample names from file list: " + fileListPath);
        TextFile filelist = new TextFile(fileListPath, false);
        sampleNames = filelist.readAsArray(0, TextFile.tab);
        filelist.close();
        System.out.println("Found " + sampleNames.length + " samples");

    }

    /**
     * Gets sample names.Takes the folder name (contained in dirName) from genotype files in the folder dirName whose names match pattern (regex)
     *
     * @param dirName
     * @param genPattern regexp pattern for matching genotype files
     *
     */
    public void loadSampleNames(String dirName, String genPattern){
        System.out.println("Getting sample names from: " + dirName + "\tfrom paths to genotype files matching: " + genPattern);
        File dir = new File(dirName);
        ArrayList<String> sNames = new ArrayList<String>();
        for (File ch : dir.listFiles())
            if (ch.isDirectory())
                for (File child : ch.listFiles())
                    if ((child.getName().matches(genPattern)) && (child.length() != 0))
                        sNames.add(ch.getName());
        sampleNames = sNames.toArray(new String[0]);
        System.out.println("Found " + sampleNames.length + " samples");
    }
}
