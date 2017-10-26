/**
 * 进度条控制器
 */
public class Scheduler {
    private int cur;
    private int max;
    private int lastPercentage;

    public Scheduler(int max){
        this.cur = 0;
        this.max = max;
        this.lastPercentage = 0;
    }

    public int getNewSchedule(){
        cur ++;
        int percentage = cur * 100/ max;
        if (percentage == lastPercentage){
            return -1;
        }else{
            lastPercentage = percentage;
            return percentage;
        }
    }
}
