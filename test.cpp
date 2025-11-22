#include <graphics.h>
#include <cmath>
#include <iostream>
#include <conio.h>
#include <Windows.h>
#include <mmsystem.h>
#include <tchar.h>
#include <stdlib.h>
#include <time.h>

struct BulletParams {
    IMAGE* imgShot;       // 图片
    int imgwidth;
    int imglength;
    int actualsize;          // 判定大小 
    double v;             // 移动速度
    double f;             // 当前帧的f值
    int centerX;          // X
    int centerY;          // Y
    double angleOffset;   // 初始角度偏移量
};

struct Character {
    IMAGE* imgChara;       // 默认贴图（朝右）
    IMAGE* imgCharaLeft;   // 左移贴图
    IMAGE* imgCharaRight;  // 右移贴图
    int imgwidth;          // 图片大小
    int imglength;
    int actualsize;
    float speed;          // 绝对速度（提高敌人速度）
    double vx;            // 速度向量(x,y)
    double vy;
    int cX;               // X,Y
    int cY; 
    bool showHitPoint;    // 是否显示判定点
    int direction;        // 方向：0=默认（朝右），1=左移，2=右移
};

//随机数种子
void InitRandomSeed() {
    static bool seeded = false;
    if (!seeded) {
        srand((unsigned int)time(NULL) + GetTickCount());
        seeded = true;
    }
}

//蝶弹自机狙 
void DrawBullet_butter(const BulletParams& params, const Character& chara, const Character& enemy, int* temp, DWORD* lastHitTime, double progress, bool isDead) {
    const int ROUND_INTERVAL = 100;         // 每轮间隔帧数
    const int ACTIVE_DURATION = 120;        // 每轮活跃时长（帧）
    const int FIRE_POINT_INTERVAL = 50;     // 发射点间距（像素）
    const int MAX_AXIS_LENGTH = 600;        // 单轴最大长度（像素）
    const int MAX_FIRE_POINTS_PER_AXIS = 20;// 单轴最大发射点数量
    const int MAX_BULLETS = 1000;           // 最大同时存在弹幕数
    const int BULLET_PER_FIRE = 3;          // 每个发射点单次发射数
    const double MAX_BULLET_DISTANCE = 1000;// 弹幕最大飞行距离
    const double PI = acos(-1.0);           
    const double BULLET_ANGLE_STEP = 2 * PI / 3;
    const double PARAM3_SPEED = MAX_AXIS_LENGTH / (double)ACTIVE_DURATION; //聚合酶移动速度

    static int fireCycle = 0;               // 计时器 
    static bool isRoundActive = false;      // 是否处于发射期
    static int activeFrameCounter = 0;      // 活跃期内帧计数
    
    //复制链属性 
    static double axisAngle1 = 0.0;        
    static double axisAngle2 = 0.0;         
    static POINT firePoints1[MAX_FIRE_POINTS_PER_AXIS]; // 轴1发射点
    static POINT firePoints2[MAX_FIRE_POINTS_PER_AXIS]; // 轴2发射点
    static int axis1PointCount = 0;         // 轴1有效发射点数量
    static int axis2PointCount = 0;         // 轴2有效发射点数量
    static bool axis1PointActivated[MAX_FIRE_POINTS_PER_AXIS] = {false}; // 轴1激活标记
    static bool axis2PointActivated[MAX_FIRE_POINTS_PER_AXIS] = {false}; // 轴2激活标记
    static bool axis1PointFired[MAX_FIRE_POINTS_PER_AXIS] = {false};     // 轴1发射标记
    static bool axis2PointFired[MAX_FIRE_POINTS_PER_AXIS] = {false};     // 轴2发射标记

	//聚合酶属性 
    static int param3_1_x = 0, param3_1_y = 0;
    static int param3_2_x = 0, param3_2_y = 0;
    static int param3Progress1 = 0, param3Progress2 = 0;
    static bool param3Active1 = false, param3Active2 = false;
    
    //普通弹幕属性 
    static int bulletFireX[MAX_BULLETS];
    static int bulletFireY[MAX_BULLETS];
    static double bulletAngle[MAX_BULLETS];
    static int bulletFireFrame[MAX_BULLETS];
    static bool bulletActive[MAX_BULLETS] = {false};
    static int bulletIndex = 0, bulletCount = 0;

    fireCycle++;
    if (fireCycle % ROUND_INTERVAL == 0) {
        isRoundActive = true;
        activeFrameCounter = 0;
        axis1PointCount = 0;
        axis2PointCount = 0;
        param3Progress1 = param3Progress2 = 0;
        param3Active1 = param3Active2 = true;

        memset(axis1PointActivated, false, sizeof(axis1PointActivated));
        memset(axis2PointActivated, false, sizeof(axis2PointActivated));
        memset(axis1PointFired, false, sizeof(axis1PointFired));
        memset(axis2PointFired, false, sizeof(axis2PointFired));

        double dx = chara.cX - enemy.cX;
        double dy = chara.cY - enemy.cY;
        double baseLineAngle = atan2(dy, dx); // 敌机->自机的连线角度
        axisAngle1 = baseLineAngle + PI / 4;  // 轴1：逆时针45°
        axisAngle2 = baseLineAngle - PI / 4;  // 轴2：顺时针45°
        // 角度处理
        if (axisAngle1 >= 2 * PI) axisAngle1 -= 2 * PI;
        if (axisAngle2 < 0) axisAngle2 += 2 * PI;

		// 聚合酶逻辑 
		for (int n = 1;; n++) {
		    double dist = n * FIRE_POINT_INTERVAL;
		    if (dist > MAX_AXIS_LENGTH || axis1PointCount >= MAX_FIRE_POINTS_PER_AXIS) break;
		    int x = (int)(enemy.cX + dist * cos(axisAngle1));
		    int y = (int)(enemy.cY + dist * sin(axisAngle1));
		    if (x >= 0 && x < 600 && y >= 0 && y < 700) { 
		        firePoints1[axis1PointCount].x = x;
		        firePoints1[axis1PointCount].y = y;
		        axis1PointCount++;
		    } else break;
		}
		for (int n = 1;; n++) {
		    double dist = n * FIRE_POINT_INTERVAL;
		    if (dist > MAX_AXIS_LENGTH || axis2PointCount >= MAX_FIRE_POINTS_PER_AXIS) break;
		    int x = (int)(enemy.cX + dist * cos(axisAngle2));
		    int y = (int)(enemy.cY + dist * sin(axisAngle2));
		    if (x >= 0 && x < 600 && y >= 0 && y < 700) {
		        firePoints2[axis2PointCount].x = x;
		        firePoints2[axis2PointCount].y = y;
		        axis2PointCount++;
		    } else break;
		}
        // 初始化聚合酶
        param3_1_x = param3_2_x = enemy.cX;
        param3_1_y = param3_2_y = enemy.cY;
    }

    // PCR逻辑 
    if (isRoundActive) {
        if (activeFrameCounter++ >= ACTIVE_DURATION) {
            isRoundActive = param3Active1 = param3Active2 = false;
            return;
        }
        if (param3Active1) {
            double dist = param3Progress1 * PARAM3_SPEED;
            param3_1_x = (int)(enemy.cX + dist * cos(axisAngle1));
            param3_1_y = (int)(enemy.cY + dist * sin(axisAngle1));
            param3Progress1++;
            if (param3_1_x < 0 || param3_1_x >= 600 || param3_1_y < 0 || param3_1_y >= 700 || dist >= MAX_AXIS_LENGTH) {
                param3Active1 = false;
            } else {
                // 激活途经发射点
                for (int i = 0; i < axis1PointCount; i++) {
                    if (!axis1PointActivated[i] && dist >= (i+1)*FIRE_POINT_INTERVAL - 1 && dist <= (i+1)*FIRE_POINT_INTERVAL + 1) {
                        axis1PointActivated[i] = true;
                    }
                }
            }
        }
        if (param3Active2) {
            double dist = param3Progress2 * PARAM3_SPEED;
            param3_2_x = (int)(enemy.cX + dist * cos(axisAngle2));
            param3_2_y = (int)(enemy.cY + dist * sin(axisAngle2));
            param3Progress2++;
            // 边界检测
            if (param3_2_x < 0 || param3_2_x >= 600 || param3_2_y < 0 || param3_2_y >= 700 || dist >= MAX_AXIS_LENGTH) {
                param3Active2 = false;
            } else {
                // 激活途经发射点
                for (int i = 0; i < axis2PointCount; i++) {
                    if (!axis2PointActivated[i] && dist >= (i+1)*FIRE_POINT_INTERVAL - 1 && dist <= (i+1)*FIRE_POINT_INTERVAL + 1) {
                        axis2PointActivated[i] = true;
                    }
                }
            }
        }


        // 绘制聚合酶
        if (param3Active1 && !isDead) {
            putimage(param3_1_x - params.imgwidth/2, param3_1_y - params.imglength/2, params.imgShot, SRCPAINT);
        }
        if (param3Active2 && !isDead) {
            putimage(param3_2_x - params.imgwidth/2, param3_2_y - params.imglength/2, params.imgShot, SRCPAINT);
        }

        // 聚合酶催化产生三条弹幕链 
        for (int i = 0; i < axis1PointCount; i++) {
            if (axis1PointActivated[i] && !axis1PointFired[i]) {
                int x = firePoints1[i].x, y = firePoints1[i].y;
                double dirAngle = atan2(chara.cY - y, chara.cX - x); // 朝向自机
                for (int j = 0; j < BULLET_PER_FIRE; j++) {
                    double angle = dirAngle + j * BULLET_ANGLE_STEP;
                    angle = angle >= 2*PI ? angle-2*PI : (angle < 0 ? angle+2*PI : angle);
                    // 存储弹幕状态
                    bulletFireX[bulletIndex] = x;
                    bulletFireY[bulletIndex] = y;
                    bulletAngle[bulletIndex] = angle;
                    bulletFireFrame[bulletIndex] = fireCycle;
                    bulletActive[bulletIndex] = true;
                    // 循环更新索引
                    bulletIndex = (bulletIndex + 1) % MAX_BULLETS;
                    if (bulletCount < MAX_BULLETS) bulletCount++;
                }
                axis1PointFired[i] = true;
            }
        }
        for (int i = 0; i < axis2PointCount; i++) {  //同上Ctrl+V 
            if (axis2PointActivated[i] && !axis2PointFired[i]) {
                int x = firePoints2[i].x, y = firePoints2[i].y;
                double dirAngle = atan2(chara.cY - y, chara.cX - x);
                for (int j = 0; j < BULLET_PER_FIRE; j++) {
                    double angle = dirAngle + j * BULLET_ANGLE_STEP;
                    angle = angle >= 2*PI ? angle-2*PI : (angle < 0 ? angle+2*PI : angle);
                    bulletFireX[bulletIndex] = x;
                    bulletFireY[bulletIndex] = y;
                    bulletAngle[bulletIndex] = angle;
                    bulletFireFrame[bulletIndex] = fireCycle;
                    bulletActive[bulletIndex] = true;
                    bulletIndex = (bulletIndex + 1) % MAX_BULLETS;
                    if (bulletCount < MAX_BULLETS) bulletCount++;
                }
                axis2PointFired[i] = true;
            }
        }
    }

    // 绘制+碰撞检测
    for (int b = 0; b < bulletCount; b++) {
        int idx = (bulletIndex - bulletCount + b + MAX_BULLETS) % MAX_BULLETS;
        if (!bulletActive[idx]) continue;
        int fireX = bulletFireX[idx], fireY = bulletFireY[idx];
        double angle = bulletAngle[idx];
        int flightFrames = fireCycle - bulletFireFrame[idx];
        double dist = params.v * flightFrames;
        if (dist > MAX_BULLET_DISTANCE) {
            bulletActive[idx] = false;
            continue;
        }
        int x = (int)(fireX + dist * cos(angle));
        int y = (int)(fireY + dist * sin(angle));

        // 绘制弹幕（屏幕内有效）
        if (x >= -params.imgwidth/2 && x <= 600 + params.imgwidth/2 && y >= -params.imglength/2 && y <= 700 + params.imglength/2) {
            if (!isDead) {
                putimage(x - params.imgwidth/2, y - params.imglength/2, params.imgShot, SRCPAINT);
            }
            // 碰撞检测
            double distToChara = sqrt(pow(x - chara.cX, 2) + pow(y - chara.cY, 2));
            if (distToChara + chara.actualsize <= params.actualsize) {
                DWORD currTime = GetTickCount();
                if (currTime - *lastHitTime >= 1000 && !isDead) {
                    (*temp)++;
                    // 播放碰撞音效
                    char cmd[512];
                    sprintf(cmd, "close hit_sound"); mciSendString(cmd, NULL, 0, NULL);
                    sprintf(cmd, "open sounds/biu.wav alias hit_sound"); mciSendString(cmd, NULL, 0, NULL);
                    sprintf(cmd, "play hit_sound"); mciSendString(cmd, NULL, 0, NULL);
                    *lastHitTime = currTime;
                }
            }
	        } else {
	            bulletActive[idx] = false;
        }
    }
}

//开花大玉弹 
void DrawBullet_yuyuko(const BulletParams& params, const Character& chara, int* temp, DWORD* lastHitTime, double progress, bool isDead) {
    const int BULLET_COUNT_PER_ROUND = 40;    // 每轮弹幕数量
    
    const int WINDOW_WIDTH = 600;
    const int WINDOW_HEIGHT = 700;
    const double ANGLE_STEP = 2 * acos(-1) / BULLET_COUNT_PER_ROUND;
    const int GROUP_INNER_INTERVAL = 15;      //组内间隔帧 
    const int FIRE_GROUP_INTERVAL = 300;      //每轮间隔帧 
    const int ROUND_COUNT_PER_GROUP = 3;     
    const double MAX_DISTANCE = WINDOW_WIDTH * 1.2; 
    const int MAX_FIRED_ROUNDS = 100;  
    static double roundBaseAngles[MAX_FIRED_ROUNDS] = {0.0};  // 基准角度
    static int roundFireFrames[MAX_FIRED_ROUNDS] = {0};       // 发射帧序号
    static int roundCount = 0;                                // 有效轮次数量
    static int roundIndex = 0;                                // 循环数组索引
    static int fireCycle = 0;                                 // 计时器 
    static bool isGroupActive = false;         
    static int innerFrameCounter = 0;          // 组内计时器
    static int roundInGroup = 0;               // 组内当前轮次(123)
    static double groupFixedAngle = 0.0;       // 实际角度
    //全局帧计数
    fireCycle++;

    if (fireCycle % FIRE_GROUP_INTERVAL == 0 && !isGroupActive) {
        isGroupActive = true;     
        innerFrameCounter = 0;     
        roundInGroup = 0;         
    }

    // 3组内发射逻辑
    if (isGroupActive) {
        innerFrameCounter++;
        if (innerFrameCounter >= GROUP_INNER_INTERVAL) {
            if (roundInGroup == 0) {
                double dx = chara.cX - params.centerX; 
                double dy = params.centerY - chara.cY;  
                groupFixedAngle = atan2(dy, dx); 
            }
            // 存储当前轮数据
            roundBaseAngles[roundIndex] = groupFixedAngle;
            roundFireFrames[roundIndex] = fireCycle;
            // 循环数组更新
            roundIndex = (roundIndex + 1) % MAX_FIRED_ROUNDS;
            if (roundCount < MAX_FIRED_ROUNDS) {
                roundCount++;
            }
            roundInGroup++;
            if (roundInGroup >= ROUND_COUNT_PER_GROUP) isGroupActive = false;
            innerFrameCounter = 0;  // 重置组内计时器 
        }
    }

    for (int r = 0; r < roundCount; r++) {
        // 计算当前轮在循环数组中的实际索引
        int currIdx = (roundIndex - roundCount + r + MAX_FIRED_ROUNDS) % MAX_FIRED_ROUNDS;
        double currAngle = roundBaseAngles[currIdx];
        int currFireFrame = roundFireFrames[currIdx];
        int flightFrames = fireCycle - currFireFrame;
        double distance = params.v * flightFrames;
        if (distance > MAX_DISTANCE) continue;
            
        // 绘制&判定当前轮的所有弹幕
        for (int i = 0; i < BULLET_COUNT_PER_ROUND; i++) {
            double bulletAngle = i * ANGLE_STEP + currAngle;
            int winX = params.centerX + static_cast<int>(distance * cos(bulletAngle));
            int winY = params.centerY - static_cast<int>(distance * sin(bulletAngle));
            if (winX >= -params.imgwidth/2 && winX <= WINDOW_WIDTH + params.imgwidth/2 &&
                winY >= -params.imglength/2 && winY <= WINDOW_HEIGHT + params.imglength/2) {
                putimage(winX - params.imgwidth/2, winY - params.imglength/2, params.imgShot, SRCPAINT);
                
                double distToChara = sqrt(pow(winX - chara.cX, 2) + pow(winY - chara.cY, 2));
                bool isCollided = (distToChara + chara.actualsize <= params.actualsize);

                DWORD currentTime = GetTickCount();
                if (isCollided && (currentTime - *lastHitTime >= 1000) && !isDead) {
                    (*temp)++;
                    LPCSTR biu_sound = "sounds/biu.wav";	
                    char cmd_close[512], cmd_open[512], cmd_play[512];
                    sprintf(cmd_close, "close hit_sound");
                    mciSendString(cmd_close, NULL, 0, NULL);
                    sprintf(cmd_open, "open %s alias hit_sound", biu_sound);
                    mciSendString(cmd_open, NULL, 0, NULL);
                    sprintf(cmd_play, "play hit_sound");
                    mciSendString(cmd_play, NULL, 0, NULL);        
                    *lastHitTime = currentTime; 
                }
            }
        }
    }
}

//波与粒的境界 
void DrawBullet_yokary(const BulletParams& params, const Character& chara, int* temp, DWORD* lastHitTime, double progress, bool isDead) {
    const int WINDOW_WIDTH = 600;
    const int WINDOW_HEIGHT = 700;
    const double MAX_DISTANCE = 1000;

    double tMax = 300 * progress;
    for (double t = 0; t <= tMax; t += 0.5) {
        double distance = params.v * t;
        if (distance > MAX_DISTANCE) continue;
        int nMax = static_cast<int>(params.f - t);
        nMax = nMax < 0 ? 0 : nMax;

        // 计算角度
        double angleSum = 0.0;
        for (int n = 0; n <= nMax; ++n) {
            double a = (acos(-1) / 180) * 0.1 * n;
            angleSum += (a + params.angleOffset);
        }

        // 计算坐标
        double x = distance * cos(angleSum);
        double y = distance * sin(angleSum);
        
        // 转换为窗口坐标
        int winX = params.centerX + static_cast<int>(x);
        int winY = params.centerY - static_cast<int>(y);

        if (!(winX >= -params.imgwidth/2 && winX <= WINDOW_WIDTH + params.imgwidth/2 &&
              winY >= -params.imglength/2 && winY <= WINDOW_HEIGHT + params.imglength/2)) continue; 

        // 绘制弹幕
        if (!isDead) putimage(winX - params.imgwidth/2, winY - params.imglength/2, params.imgShot, SRCPAINT);

        // 碰撞检测
        double distToChara = sqrt(pow(winX - chara.cX, 2) + pow(winY - chara.cY, 2));
        bool isCollided = (distToChara + chara.actualsize <= params.actualsize);

        // 无敌帧判断
        DWORD currentTime = GetTickCount();
        if (isCollided && (currentTime - *lastHitTime >= 1000) && !isDead) {
            (*temp)++;         
            LPCSTR biu_sound = "sounds/biu.wav";	
            char cmd_close[512];
            char cmd_open[512];
            char cmd_play[512];
            sprintf(cmd_close, "close hit_sound");
            mciSendString(cmd_close, NULL, 0, NULL);
            sprintf(cmd_open, "open %s alias hit_sound", biu_sound);
            mciSendString(cmd_open, NULL, 0, NULL);
            sprintf(cmd_play, "play hit_sound");
            mciSendString(cmd_play, NULL, 0, NULL);        
            *lastHitTime = currentTime; 
        }
    }
}

// 灰度效果
void ApplyGrayEffect(IMAGE* img) {
    DWORD* pBuf = GetImageBuffer(img);
    int width = img->getwidth();
    int height = img->getheight();
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            DWORD color = pBuf[y * width + x];
            BYTE r = GetRValue(color);
            BYTE g = GetGValue(color);
            BYTE b = GetBValue(color);
            BYTE gray = (BYTE)(0.299 * r + 0.587 * g + 0.114 * b);
            pBuf[y * width + x] = RGB(gray, gray, gray);
        }
    }
}

// 敌机移动逻辑
void UpdateEnemyMovement(Character& enemy) {
    static int frameCounter = 0;          // 计时器 
    static int targetX = 0, targetY = 0;  // 目标位置
    static double moveDirX = 0, moveDirY = 0; // 移动方向向量
    frameCounter++;
    if (frameCounter % 500 == 0) {
        targetX = 100 + rand() % (350);
        targetY = 75 + rand() % (25);

        double dx = targetX - enemy.cX;
        double dy = targetY - enemy.cY;
        double dist = sqrt(dx * dx + dy * dy);
        if (dist > 0) {
            moveDirX = dx / dist;
            moveDirY = dy / dist;
        }
        frameCounter = 0;
    }

    int newX = enemy.cX + (int)(enemy.speed * moveDirX);
    int newY = enemy.cY + (int)(enemy.speed * moveDirY);
    if (newX >= 100 && newX <= 500) enemy.cX = newX;
    if (newY >= 75 && newY <= 150) enemy.cY = newY;

    if (moveDirX < -0.5) enemy.direction = 1;   
    else if (moveDirX > 0.5) enemy.direction = 2;
    if (enemy.actualsize == enemy.cX) enemy.direction = 0;
}

int main() {
    // 常量设置
    const double F_MAX = 5000.0;
    const double PI = acos(-1);
    const int SLEEP_MS = 10;
    const int WINDOW_WIDTH = 600;
    const int WINDOW_HEIGHT = 700;
    ExMessage msg = { 0 };
    int temp = 0; 
    DWORD lastHitTime = 0;
    DWORD startTime = 0;
    DWORD startBulletTime = 0;
    const int FADE_DURATION = 3500; 
    bool isDead = false;
    bool deadMusicPlayed = false; 
  
    // 初始化随机数种子
    InitRandomSeed();

    // 初始化窗口
    IMAGE background,background1,leftbackground,lives;
    loadimage(&background, "assest//eff06b.png", WINDOW_WIDTH, WINDOW_HEIGHT);
    loadimage(&background1, "assest//eff06.png", WINDOW_WIDTH, WINDOW_HEIGHT);
	loadimage(&leftbackground, "assest//leftbackground.png", 300, WINDOW_HEIGHT);
	loadimage(&lives, "assest//lives.png", 30, 30); 
    
    initgraph(WINDOW_WIDTH+300, WINDOW_HEIGHT, EX_DBLCLKS | EW_SHOWCONSOLE);
    cleardevice();
    startTime = GetTickCount();

    // 弹幕角度设置
    const double angleOffsets[] = {8 * PI / 5, 6 * PI /5, 4 * PI / 5, 2 * PI / 5, 0.0};
    const int BULLET_GROUP_COUNT = sizeof(angleOffsets) / sizeof(angleOffsets[0]);
    
    // 初始化弹幕参数
    IMAGE Param1Shot,Param2Shot,Param3Shot;
    BulletParams param1 = {
        &Param1Shot,    // imgShot
        8,          // Paramwidth
        16,          // paramlength
        8,          // paramsize(actual)
        2.5,         // v
        100.0,       // f
        WINDOW_WIDTH / 2,    // centerX
        WINDOW_HEIGHT / 2,    // centerY
        0.0          // angle
    };//米弹 
    BulletParams param2 = {
        &Param2Shot,    // imgShot
        50,          // Paramwidth
        50,          // paramlength
        20,          // paramsize(actual)
        1,         // v
        100.0,       // f
        WINDOW_WIDTH / 2,    // centerX
        WINDOW_HEIGHT / 2,    // centerY
        0.0          // angle
    };//大玉弹 
    BulletParams param3 = {
        &Param3Shot,    // imgShot
        30,          // Paramwidth
        30,          // paramlength
        12,          // paramsize(actual)
        2,         // v
        100.0,       // f
        WINDOW_WIDTH / 2,    // centerX
        WINDOW_HEIGHT / 2,    // centerY
        0.0          // angle
    };//蝶弹 

    loadimage(&Param1Shot, "assest//smallshot.png", param1.imgwidth, param1.imglength);
    loadimage(&Param2Shot, "assest//bigshot.png", param2.imgwidth, param2.imglength);
    loadimage(&Param3Shot, "assest//buttershot.png", param3.imgwidth, param3.imglength);

    // 初始化角色和敌人贴图
    IMAGE imgChara,imgCharaLeft,imgCharaRight;
    IMAGE imgEnemy, imgEnemyLeft, imgEnemyRight, dead;
    
    Character chara1 = {
        &imgChara,           // imgChara
        &imgCharaLeft,       // imgCharaLeft
        &imgCharaRight,      // imgCharaRight
        48,                  // imgwidth
        32,                  // imglength
        3,                   // actualsize
        2,                   // speed
        0,                   // vx
        0,                   // vy
        WINDOW_WIDTH / 2,    // cX
        WINDOW_HEIGHT - 32,  // cY
        false,               // showHitPoint
        0                    // direction
    };
    
    Character enemy1 ={
		&imgEnemy,           // 默认贴图
	    &imgEnemyLeft,       // 左移贴图
	    &imgEnemyRight,      // 右移贴图
	    144,                  // imgwidth
	    84,                  // imglength
	    0,                   // actualsize
	    3,                  
	    0,                   // vx
	    0,                   // vy
	    300,                  // 初始X坐标
	    150,                  // 初始Y坐标
	    false,               // showHitPoint
	    0                  
    };
    
    // 加载角色和敌人图片
    loadimage(&imgChara, "assest//player.png", chara1.imglength, chara1.imgwidth);
    loadimage(&imgCharaLeft, "assest//playerleft.png", chara1.imglength, chara1.imgwidth);
    loadimage(&imgCharaRight, "assest//playerright.png", chara1.imglength, chara1.imgwidth);
    loadimage(&imgEnemy, "assest//enemy.png", enemy1.imglength, enemy1.imgwidth); 
    loadimage(&imgEnemyLeft, "assest//enemyleft.png", enemy1.imglength, enemy1.imgwidth);
    loadimage(&imgEnemyRight, "assest//enemyright.png", enemy1.imglength, enemy1.imgwidth);
    loadimage(&dead, "assest//dead.png", enemy1.imglength, enemy1.imgwidth);
    
    // 播放背景音乐
    LPCSTR music_name = "sounds//background_music.mp3";
    char cmd_open[512];
    sprintf(cmd_open, "open %s alias bgm", music_name);
    mciSendString(cmd_open, NULL, 0, NULL);
    mciSendString("play bgm", NULL, 0, NULL);
    
    // 创建灰度背景副本
    IMAGE grayBackground, grayBackground1;
    grayBackground = background;
    grayBackground1 = background1;
    
    // 主循环
    while (param1.f <= F_MAX && temp <= 5) {
        BeginBatchDraw();  
        cleardevice();
		
		// 角色是否死亡 
        if (temp >= 5 && !isDead) { 
            isDead = true;
            // 应用灰度效果
            ApplyGrayEffect(&grayBackground);
            ApplyGrayEffect(&grayBackground1);
            
            // 播放死亡音效 
            if (!deadMusicPlayed) {
                mciSendString("close bgm", NULL, 0, NULL);
                LPCSTR dead_music_name = "sounds//dead.mp3";
                char deadmusic_open[512];
                sprintf(deadmusic_open, "open %s alias deadbgm", dead_music_name);
                mciSendString(deadmusic_open, NULL, 0, NULL);
                mciSendString("play deadbgm", NULL, 0, NULL);
                deadMusicPlayed = true;
            }
        }
        
        // 绘制背景
        if (isDead) {
            putimage(0, 0, &grayBackground);
            putimage(0, 0, &grayBackground1, SRCCOPY);
        } else {
            putimage(0, 0, &background);
            putimage(0, 0, &background1, SRCPAINT);
        }

        // 处理键盘消息
        ExMessage msg = {0};
        while (peekmessage(&msg, EX_KEY)) {
            if (msg.message == WM_KEYDOWN) {
                switch (msg.vkcode) {
                    case VK_UP:    chara1.vy = -1.5; break;
                    case VK_DOWN:  chara1.vy = 2;  break;
                    case VK_LEFT:  
                        chara1.vx = -1.5; 
                        chara1.direction = 1;
                        break;
                    case VK_RIGHT: 
                        chara1.vx = 1.75;  
                        chara1.direction = 2; 
                        break;
                    case VK_SHIFT: 
                        chara1.speed = 1; 
                        chara1.showHitPoint = true;
                        break;
                }
            } else if (msg.message == WM_KEYUP) {
                switch (msg.vkcode) {
                    case VK_UP:    chara1.vy = 0; break;
                    case VK_DOWN:  chara1.vy = 0; break;
                    case VK_LEFT:  
                        chara1.vx = 0; 
                        if (chara1.direction == 1) chara1.direction = 0;
                        break;
                    case VK_RIGHT: 
                        chara1.vx = 0; 
                        if (chara1.direction == 2) chara1.direction = 0;
                        break;
                    case VK_SHIFT: 
                        chara1.speed = 2; 
                        chara1.showHitPoint = false;
                        break;
                }
            }
        }

        // 角色移动 + 边界检测
        chara1.cX += chara1.speed * chara1.vx;
        chara1.cY += chara1.speed * chara1.vy;  
        if (chara1.cX - chara1.imglength/2 < 0) chara1.cX = chara1.imglength/2;
        if (chara1.cX + chara1.imglength/2 > WINDOW_WIDTH) chara1.cX = WINDOW_WIDTH - chara1.imglength/2;
        if (chara1.cY - chara1.imgwidth/2 < 0) chara1.cY = chara1.imgwidth/2;
        if (chara1.cY + chara1.imgwidth/2 > WINDOW_HEIGHT) chara1.cY = WINDOW_HEIGHT - chara1.imgwidth/2;
        enemy1.actualsize = enemy1.cX;
        if (!isDead) UpdateEnemyMovement(enemy1);

        // 弹幕绘制逻辑
        DWORD currentTime = GetTickCount();
        if (currentTime - startTime >= 5000 && !isDead) {
            if (startBulletTime == 0) startBulletTime = currentTime;

            // 时间判断 
            double progress = static_cast<double>(currentTime - startBulletTime) / FADE_DURATION;
            progress = progress > 1.0 ? 1.0 : progress; 
            progress = progress < 0.0 ? 0.0 : progress; 

            // 绘制弹幕
            for (int i = 0; i < BULLET_GROUP_COUNT; ++i) {
                param1.angleOffset = angleOffsets[i];
                DrawBullet_yokary(param1, chara1, &temp, &lastHitTime, progress, isDead);
            }
            DrawBullet_yuyuko(param2, chara1, &temp, &lastHitTime, progress, isDead);
            DrawBullet_butter(param3, chara1, enemy1, &temp, &lastHitTime, progress, isDead);
        }

        // 绘制机体 
        if (!isDead){
	        IMAGE* enemyCurrentImg = enemy1.imgChara; // 默认贴图
	        if (enemy1.direction == 1) {
	            enemyCurrentImg = enemy1.imgCharaLeft; // 左移贴图
	        } else if (enemy1.direction == 2) {
	            enemyCurrentImg = enemy1.imgCharaRight; // 右移贴图
	        }
	        putimage(enemy1.cX - enemy1.imglength/2, enemy1.cY - enemy1.imgwidth/2, enemyCurrentImg, SRCPAINT);
		}
        if (!isDead && (currentTime - lastHitTime >= 1000 || (currentTime % 200) < 100)) {
            IMAGE* currentImg = chara1.imgChara;
            if (chara1.direction == 1) { 
                currentImg = chara1.imgCharaLeft;
            } else if (chara1.direction == 2) { 
                currentImg = chara1.imgCharaRight;
            }
            putimage(chara1.cX - chara1.imglength/2, chara1.cY - chara1.imgwidth/2, currentImg, SRCPAINT);
            if (chara1.showHitPoint) {
                setfillcolor(WHITE);
                setlinecolor(RED);
                fillcircle(chara1.cX, chara1.cY, chara1.actualsize);
            }
        }
        
        // 绘制生命数
        putimage(600, 0, &leftbackground, SRCCOPY);
		for (int i=0; i<5-temp; i++){
			putimage(700+i*30, 123, &lives);
		}
        
        // 显示死亡信息
        if (isDead) {
            settextcolor(WHITE);
            setbkmode(TRANSPARENT);
            outtextxy(WINDOW_WIDTH/2 - 50, WINDOW_HEIGHT/2, _T("满  身  疮  痍"));
            putimage(WINDOW_WIDTH/2 - 50, 100, &dead, SRCPAINT);
        }
        
        EndBatchDraw(); 
        Sleep(SLEEP_MS);
        param1.f += 1.0;
    }
    if (!isDead) mciSendString("close bgm", NULL, 0, NULL);

    _getch();
    closegraph();
    return 0;
}
