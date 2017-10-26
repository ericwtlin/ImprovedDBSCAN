package utils;

import java.io.File;

public class FileOperation {

    public static void mkDir(String dirPath){
        File dir = new File(dirPath);
        if (!dir.exists()){
            if (dir.mkdirs()){
                System.out.println(String.format("Dir %s created", dirPath));
            }else{
                System.out.println(String.format("Create dir %s failed", dirPath));
            }
        }else{
            //throw new Exception(String.format("Dir %s already exists!", groupDirPath));
            FileOperation.cleanDir(dir);
        }
        try {
            Thread.sleep(3000);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }
	public static boolean cleanDir(File dir){
       if (dir.isDirectory()) {
            String[] children = dir.list();
            for (int i=0; i<children.length; i++) {
                boolean success = deleteDir(new File(dir, children[i]));
                if (!success) {
                    return false;
                }
            }
        }
        return true;
    }
	public static boolean deleteDir(File dir) {
        if (dir.isDirectory()) {
            String[] children = dir.list();
            for (int i=0; i<children.length; i++) {
                boolean success = deleteDir(new File(dir, children[i]));
                if (!success) {
                    return false;
                }
            }
        }
        // 目录此时为空，可以删除
        return dir.delete();
    }
}
