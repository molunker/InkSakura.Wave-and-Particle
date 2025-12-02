#include <graphics.h>
#include <cmath>
#include <iostream>
#include <conio.h>
#include <Windows.h>
#include <mmsystem.h>
#include <vector>
#include <numeric>

#pragma comment(lib, "winmm.lib")

const int WINDOW_WIDTH = 600;
const int WINDOW_HEIGHT = 700;
const double PI = acos(-1.0);

struct BulletParams {
    IMAGE* imgShot;
    int imgwidth;
    int imglength;
    int actualsize;
    double v;
    double f;
    int centerX;
    int centerY;
    double angleOffset;
};

struct Character {
    IMAGE* imgChara;
    IMAGE* imgCharaLeft;
    IMAGE* imgCharaRight;
    int imgwidth;
    int imglength;
    int actualsize;
    float speed;
    double vx;
    double vy;
    int cX;
    int cY;
    bool showHitPoint;
    int direction; // 0=default, 1=left, 2=right
};

namespace GameStates {
    // --- DrawScatterBullet 状态 ---
    struct ScatterState {
        double frame;
    };
    void ResetScatterState(ScatterState& state) {
        state.frame = 0.0;
    }

    // --- DrawBullet_butter 状态 ---
    struct ButterState {
        int fireCycle;
        bool isRoundActive;
        int activeFrameCounter;
        double axisAngle1;
        double axisAngle2;
        static const int MAX_FIRE_POINTS_PER_AXIS = 20;
        POINT firePoints1[MAX_FIRE_POINTS_PER_AXIS];
        POINT firePoints2[MAX_FIRE_POINTS_PER_AXIS];
        int axis1PointCount;
        int axis2PointCount;
        bool axis1PointActivated[MAX_FIRE_POINTS_PER_AXIS];
        bool axis2PointActivated[MAX_FIRE_POINTS_PER_AXIS];
        bool axis1PointFired[MAX_FIRE_POINTS_PER_AXIS];
        bool axis2PointFired[MAX_FIRE_POINTS_PER_AXIS];
        int param3_1_x, param3_1_y;
        int param3_2_x, param3_2_y;
        int param3Progress1, param3Progress2;
        bool param3Active1, param3Active2;
        static const int MAX_BULLETS = 300;
        int bulletFireX[MAX_BULLETS];
        int bulletFireY[MAX_BULLETS];
        double bulletAngle[MAX_BULLETS];
        int bulletFireFrame[MAX_BULLETS];
        bool bulletActive[MAX_BULLETS];
        int bulletIndex, bulletCount;
    };
    void ResetButterState(ButterState& state) {
        memset(&state, 0, sizeof(ButterState)); 
    }

    // --- DrawBullet_yuyuko 状态 ---
    struct YuyukoState {
        static const int MAX_FIRED_ROUNDS = 50;
        double roundBaseAngles[MAX_FIRED_ROUNDS];
        int roundFireFrames[MAX_FIRED_ROUNDS];
        int roundCount;
        int roundIndex;
        int fireCycle;
        bool isGroupActive;
        int innerFrameCounter;
        int roundInGroup;
        double groupFixedAngle;
    };
    void ResetYuyukoState(YuyukoState& state) {
        memset(&state, 0, sizeof(YuyukoState));
    }
    struct TesseractState {
        struct Vec4 { double x, y, z, w; };
        struct Vec3 { double x, y, z; };
        struct Vec2 { double x, y; };

        double theta, phi, psi;
        Vec2 projectedVertices[16];

        const Vec4 vertices[16] = {
            {1,1,1,1}, {-1,1,1,1}, {1,-1,1,1}, {-1,-1,1,1},
            {1,1,-1,1}, {-1,1,-1,1}, {1,-1,-1,1}, {-1,-1,-1,1},
            {1,1,1,-1}, {-1,1,1,-1}, {1,-1,1,-1}, {-1,-1,1,-1},
            {1,1,-1,-1}, {-1,1,-1,-1}, {1,-1,-1,-1}, {-1,-1,-1,-1}
        };
    };
    void ResetTesseractState(TesseractState& state) {
        state.theta = state.phi = state.psi = 0.0;
    }
}

// 散乱银河 (欧拉函数)
void DrawScatterBullet(GameStates::ScatterState& state, const BulletParams& params, const Character& chara, int* temp, DWORD* lastHitTime, double progress, bool isDead) {
    if (!isDead) {
         state.frame++;
    }

    const double BASE_LOG_SCALE = 0.4;
    const double BREATH_AMP = 0.35;
    const double DRAW_SCALE = 22.0; 
    const double ROTATION_SPEED = 0.005;
    const double P1_BREATH_SPEED = 0.013;
    const double P2_BREATH_SPEED = 0.017;

    double global_theta = state.frame * ROTATION_SPEED;
    double r1 = BASE_LOG_SCALE + sin(state.frame * P1_BREATH_SPEED) * BREATH_AMP;
    double r2 = BASE_LOG_SCALE + sin(state.frame * P2_BREATH_SPEED) * BREATH_AMP;
    double scale1 = exp(r1);
    double scale2 = exp(r2);
    double theta1 = 0.0 + global_theta;
    double theta2 = (PI / 2.0) + global_theta;
    double p1_x = scale1 * cos(theta1);
    double p1_y = scale1 * sin(theta1);
    double p2_x = scale2 * cos(theta2);
    double p2_y = scale2 * sin(theta2);

    const int RANGE = 15;
    for (int i = -RANGE; i <= RANGE; i++) {
        for (int j = -RANGE; j <= RANGE; j++) {
            double X_world = i * p1_x + j * p2_x;
            double Y_world = i * p1_y + j * p2_y;
            int winX = (int)(X_world * DRAW_SCALE) + WINDOW_WIDTH / 2;
            int winY = (int)(-Y_world * DRAW_SCALE) + WINDOW_HEIGHT / 2;

            if (winX < 0 || winX > WINDOW_WIDTH || winY < 0 || winY > WINDOW_HEIGHT) continue;
            if (!isDead) putimage(winX - params.imgwidth / 2, winY - params.imglength / 2, params.imgShot, SRCPAINT);

            double distToChara = sqrt(pow(winX - chara.cX, 2) + pow(winY - chara.cY, 2));
            if (distToChara < (chara.actualsize + params.actualsize) / 2.0) {
                DWORD currTime = GetTickCount();
                if (currTime - *lastHitTime >= 1514 && !isDead) {
                    (*temp)++;
                    mciSendString("play sounds/biu.wav from 0", NULL, 0, NULL);
                    *lastHitTime = currTime;
                }
            }
        }
    }
}

// 碟弹自机狙 
void DrawBullet_butter(GameStates::ButterState& state, const BulletParams& params, const Character& chara, const Character& enemy, int* temp, DWORD* lastHitTime, bool isDead) {
    const int ROUND_INTERVAL = 240; 
    const int ACTIVE_DURATION = 120;
    const int FIRE_POINT_INTERVAL = 50;
    const int MAX_AXIS_LENGTH = 1200;
    const int BULLET_PER_FIRE = 3;
    const double BULLET_ANGLE_STEP = 2 * PI / 3;
    const double PARAM3_SPEED = ((double)MAX_AXIS_LENGTH / ACTIVE_DURATION) / 3.0;

    state.fireCycle++;

    if (state.fireCycle % ROUND_INTERVAL == 1) {
        state.isRoundActive = true;
        state.activeFrameCounter = 0;
        state.axis1PointCount = 0;
        state.axis2PointCount = 0;
        state.param3Progress1 = 0;
        state.param3Progress2 = 0;
        state.param3Active1 = true;
        state.param3Active2 = true;
        memset(state.axis1PointActivated, false, sizeof(state.axis1PointActivated));
        memset(state.axis2PointActivated, false, sizeof(state.axis2PointActivated));
        memset(state.axis1PointFired, false, sizeof(state.axis1PointFired));
        memset(state.axis2PointFired, false, sizeof(state.axis2PointFired));

        double dx = chara.cX - enemy.cX;
        double dy = chara.cY - enemy.cY;
        double baseLineAngle = atan2(dy, dx);
        state.axisAngle1 = baseLineAngle + PI / 4.0;
        state.axisAngle2 = baseLineAngle - PI / 4.0;
		
		//聚合酶结合点位计算 
        auto generateFirePoints = [&](POINT* points, int& count, double angle) {
            for (int n = 1; ; n++) {
                double d = n * FIRE_POINT_INTERVAL;
                if (d > MAX_AXIS_LENGTH || count >= state.MAX_FIRE_POINTS_PER_AXIS) break;
                int x = (int)(enemy.cX + d * cos(angle));
                int y = (int)(enemy.cY + d * sin(angle));
                if (x >= 0 && x < WINDOW_WIDTH && y >= 0 && y < WINDOW_HEIGHT) {
                    points[count++] = { x, y };
                } else break;
            }
        };
        generateFirePoints(state.firePoints1, state.axis1PointCount, state.axisAngle1);
        generateFirePoints(state.firePoints2, state.axis2PointCount, state.axisAngle2);
    }

	//PCR逻辑 
    if (state.isRoundActive) {
        if (state.activeFrameCounter++ >= ACTIVE_DURATION * 3.0) { 
            state.isRoundActive = false;
        }

        //更新聚合酶位置 
        if (state.param3Active1) {
            double dist = ++state.param3Progress1 * PARAM3_SPEED;
            state.param3_1_x = (int)(enemy.cX + dist * cos(state.axisAngle1));
            state.param3_1_y = (int)(enemy.cY + dist * sin(state.axisAngle1));
            if (dist > MAX_AXIS_LENGTH){
                state.param3Active1 = false;
            } else {
                 if(!isDead) putimage(state.param3_1_x - params.imgwidth/2, state.param3_1_y - params.imglength/2, params.imgShot, SRCINVERT);
                for (int i = 0; i < state.axis1PointCount; i++) {
                    if (!state.axis1PointActivated[i] && dist >= (i + 1) * FIRE_POINT_INTERVAL * 0.9) { // 留出一点余量
                        state.axis1PointActivated[i] = true;
                    }
                }
            }
        }
        if (state.param3Active2) {
            double dist = ++state.param3Progress2 * PARAM3_SPEED;
            state.param3_2_x = (int)(enemy.cX + dist * cos(state.axisAngle2));
            state.param3_2_y = (int)(enemy.cY + dist * sin(state.axisAngle2));
            if (dist > MAX_AXIS_LENGTH){
                state.param3Active2 = false;
            } else {
                if(!isDead) putimage(state.param3_2_x - params.imgwidth/2, state.param3_2_y - params.imglength/2, params.imgShot, SRCINVERT);
                for (int i = 0; i < state.axis2PointCount; i++) {
                    if (!state.axis2PointActivated[i] && dist >= (i + 1) * FIRE_POINT_INTERVAL * 0.9) {
                        state.axis2PointActivated[i] = true;
                    }
                }
            }
        }

        // DNA延长生成 
        auto fireBullets = [&](bool* activated, bool* fired, int count, POINT* points) {
            for(int i=0; i<count; i++) if(activated[i] && !fired[i]){
                double dirAngle = atan2(chara.cY - points[i].y, chara.cX - points[i].x);
                for(int j=0; j<BULLET_PER_FIRE; j++) {
                    state.bulletFireX[state.bulletIndex] = points[i].x; state.bulletFireY[state.bulletIndex] = points[i].y;
                    state.bulletAngle[state.bulletIndex] = dirAngle + j * BULLET_ANGLE_STEP;
                    state.bulletFireFrame[state.bulletIndex] = state.fireCycle;
                    state.bulletActive[state.bulletIndex] = true;
                    state.bulletIndex = (state.bulletIndex + 1) % state.MAX_BULLETS;
                    if(state.bulletCount < state.MAX_BULLETS) state.bulletCount++;
                }
                fired[i]=true;
            }
        };
        fireBullets(state.axis1PointActivated, state.axis1PointFired, state.axis1PointCount, state.firePoints1);
        fireBullets(state.axis2PointActivated, state.axis2PointFired, state.axis2PointCount, state.firePoints2);
    }

    for (int i = 0; i < state.bulletCount; i++) {
        if (!state.bulletActive[i]) continue;
        double dist = params.v * (state.fireCycle - state.bulletFireFrame[i]);
        if (dist > 2500) { state.bulletActive[i] = false; continue; }

        int x = (int)(state.bulletFireX[i] + dist * cos(state.bulletAngle[i]));
        int y = (int)(state.bulletFireY[i] + dist * sin(state.bulletAngle[i]));

		//碰撞检测 
        if (x > -50 && x < WINDOW_WIDTH + 50 && y > -50 && y < WINDOW_HEIGHT + 50) {
            if (!isDead) putimage(x - params.imgwidth / 2, y - params.imglength / 2, params.imgShot, SRCPAINT);
            double distToChara = sqrt(pow(x - chara.cX, 2) + pow(y - chara.cY, 2));
            if (distToChara < (chara.actualsize + params.actualsize) / 2.0 ) {
                DWORD currTime = GetTickCount();
                if (currTime - *lastHitTime >= 1514 && !isDead) {
                    (*temp)++; 
                    mciSendStringA("play sounds/biu.wav from 0", NULL, 0, NULL);
                    *lastHitTime = currTime;
                }
            }
        } else {
            state.bulletActive[i] = false;//回收资源 
        }
    }
}


// 大玉开花伪随机弹 
void DrawBullet_yuyuko(GameStates::YuyukoState& state, const BulletParams& params, const Character& chara, int* temp, DWORD* lastHitTime, double progress, bool isDead) {
    const int BULLET_COUNT_PER_ROUND = 40;
    const double ANGLE_STEP = 2 * PI / BULLET_COUNT_PER_ROUND;
    const int GROUP_INNER_INTERVAL = 25;
    const int FIRE_GROUP_INTERVAL = 500;
    const int ROUND_COUNT_PER_GROUP = 3;
    const double MAX_DISTANCE = WINDOW_WIDTH * 1.5;

    state.fireCycle++;
    if (state.fireCycle % FIRE_GROUP_INTERVAL == 0 && !state.isGroupActive) {
        state.isGroupActive = true;
        state.innerFrameCounter = 0;
        state.roundInGroup = 0;
    }
	
	//第一颗大玉永远指向自机 
    if (state.isGroupActive) {
        if (++state.innerFrameCounter >= GROUP_INNER_INTERVAL) {
            if (state.roundInGroup == 0) { 
                state.groupFixedAngle = atan2(params.centerY - chara.cY, chara.cX - params.centerX);
            }
            state.roundBaseAngles[state.roundIndex] = state.groupFixedAngle;
            state.roundFireFrames[state.roundIndex] = state.fireCycle;
            state.roundIndex = (state.roundIndex + 1) % GameStates::YuyukoState::MAX_FIRED_ROUNDS;
            if (state.roundCount < GameStates::YuyukoState::MAX_FIRED_ROUNDS) state.roundCount++;

            if (++state.roundInGroup >= ROUND_COUNT_PER_GROUP) state.isGroupActive = false;
            state.innerFrameCounter = 0;
        }
    }

    // 绘制所有激活的弹幕轮
    for (int r = 0; r < state.roundCount; r++) {
        int currIdx = (state.roundIndex - state.roundCount + r + GameStates::YuyukoState::MAX_FIRED_ROUNDS) % GameStates::YuyukoState::MAX_FIRED_ROUNDS;
        int flightFrames = state.fireCycle - state.roundFireFrames[currIdx];
        double distance = params.v * flightFrames;
        if (distance > MAX_DISTANCE) continue;

        double currBaseAngle = state.roundBaseAngles[currIdx];
        for (int i = 0; i < BULLET_COUNT_PER_ROUND; i++) {
            double bulletAngle = i * ANGLE_STEP + currBaseAngle;
            int winX = params.centerX + static_cast<int>(distance * cos(bulletAngle));
            int winY = params.centerY - static_cast<int>(distance * sin(bulletAngle));

            if (winX > -20 && winX < WINDOW_WIDTH+20 && winY > -20 && winY < WINDOW_HEIGHT + 20) {
                if (!isDead) putimage(winX - params.imgwidth/2, winY - params.imglength/2, params.imgShot, SRCPAINT);
                double distToChara = sqrt(pow(winX - chara.cX, 2) + pow(winY - chara.cY, 2));
                 if (distToChara < (chara.actualsize + params.actualsize) / 2.0 ) {
                    DWORD currentTime = GetTickCount();
                    if (currentTime - *lastHitTime >= 1514 && !isDead) { 
                        (*temp)++;
                        mciSendString("play sounds/biu.wav from 0", NULL, 0, NULL);
                        *lastHitTime = currentTime;
                    }
                }
            }
        }
    }
}

// 波与粒的境界
void DrawBullet_yokary(const BulletParams& params, const Character& chara, int* temp, DWORD* lastHitTime, double progress, bool isDead, int TotalRound) {
    if(progress<=0) return;
    const double MAX_DISTANCE = 1200;
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

        if (winX > -20 && winX < WINDOW_WIDTH+20 && winY > -20 && winY < WINDOW_HEIGHT + 20)
        {
            if (!isDead) putimage(winX - params.imgwidth/2, winY - params.imglength/2, params.imgShot, SRCPAINT);
            double distToChara = sqrt(pow(winX - chara.cX, 2) + pow(winY - chara.cY, 2));

            if (distToChara < (chara.actualsize + params.actualsize) / 2.0) {
                DWORD currentTime = GetTickCount();
                if (currentTime - *lastHitTime >= 1514 && !isDead) {
                    (*temp)++;
                    mciSendString("play sounds/biu.wav from 0", NULL, 0, NULL);
                    *lastHitTime = currentTime;
                }
            }
        }
    }
}

// 4D超立方体
void DrawDemo(GameStates::TesseractState& state, const BulletParams& params, const Character& chara, int* temp, DWORD* lastHitTime, bool isDead) {
    const double ROT_SPEED = 0.0035;       // 旋转速度
    const double BULLET_INTERVAL = 35.0; // 边上的弹幕间距

    state.theta += ROT_SPEED;
    state.phi += ROT_SPEED * 0.8;
    state.psi += ROT_SPEED * 0.6;

    const int tesseractEdges[32][2] = {
        {0,1},{0,2},{0,4},{1,3},{1,5},{2,3},{2,6},{3,7},{4,5},{4,6},{5,7},{6,7},
        {8,9},{8,10},{8,12},{9,11},{9,13},{10,11},{10,14},{11,15},{12,13},{12,14},{13,15},{14,15},
        {0,8},{1,9},{2,10},{3,11},{4,12},{5,13},{6,14},{7,15}
    };

    // 渲染逻辑
    for (int i = 0; i < 16; ++i) {
        GameStates::TesseractState::Vec4 v = state.vertices[i];

        // 4D旋转
        double cosT=cos(state.theta), sinT=sin(state.theta); v = {v.x*cosT-v.y*sinT, v.x*sinT+v.y*cosT, v.z, v.w};
        double cosP=cos(state.phi),  sinP=sin(state.phi);  v = {v.x*cosP-v.z*sinP, v.y, v.x*sinP+v.z*cosP, v.w};
        double cosS=cos(state.psi),  sinS=sin(state.psi);  v = {v.x*cosS-v.w*sinS, v.y, v.z, v.x*sinS+v.w*cosS};

        // 4D到3D投影
        double p_w = v.w + 5.0; // 投影距离
        GameStates::TesseractState::Vec3 v3 = {v.x * 3.0 / p_w, v.y * 3.0 / p_w, v.z * 3.0 / p_w};

        double screenDist = 8.0;
        double scale = screenDist / (8.0 - v3.z);
        double final_draw_scale = 280.0;

        state.projectedVertices[i] = {
            v3.x * scale * final_draw_scale + params.centerX,
            -v3.y * scale * final_draw_scale + params.centerY
        };
    }

    // 绘制边和弹幕
    for (int e = 0; e < 32; ++e) {
        GameStates::TesseractState::Vec2 p1 = state.projectedVertices[tesseractEdges[e][0]];
        GameStates::TesseractState::Vec2 p2 = state.projectedVertices[tesseractEdges[e][1]];
        double dx = p2.x - p1.x;
        double dy = p2.y - p1.y;
        double edgeLength = sqrt(dx*dx + dy*dy);
        if (edgeLength < 1.0) continue;

        int bulletCount = static_cast<int>(edgeLength / BULLET_INTERVAL);
        if (bulletCount < 1) continue;

        for (int i = 0; i <= bulletCount; ++i) {
            double k = (double)i / bulletCount;
            double x = p1.x + k * dx;
            double y = p1.y + k * dy;

            if (!isDead && (x - params.imgwidth / 2.0)<600) {
                putimage(static_cast<int>(x - params.imgwidth / 2.0),
                         static_cast<int>(y - params.imglength / 2.0),
                         params.imgShot, SRCPAINT);
            }

            double distToChara = sqrt(pow(x - chara.cX, 2) + pow(y - chara.cY, 2));
            if (distToChara < (chara.actualsize + params.actualsize) / 2.0) {
                DWORD currTime = GetTickCount();
                if (currTime - *lastHitTime >= 1514 && !isDead) {
                    (*temp)++;
                    mciSendString("play sounds/biu.wav from 0", NULL, 0, NULL);
                    *lastHitTime = currTime;
                }
            }
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

void DrawSingleAsciiChar(int destX, int destY, char c, int size, IMAGE& g_AsciiFontImg) {
    const int ORIGIN_CHAR_SIZE = 16;
    unsigned char uc = static_cast<unsigned char>(c) - 32;
    int col = uc % 16;
    int row = uc / 16;
    int fontX = col * ORIGIN_CHAR_SIZE;
    int fontY = row * ORIGIN_CHAR_SIZE;
    IMAGE tempCharImg;
    SetWorkingImage(&g_AsciiFontImg);
    getimage(&tempCharImg, fontX, fontY, ORIGIN_CHAR_SIZE, ORIGIN_CHAR_SIZE);
    SetWorkingImage();
    IMAGE scaledCharImg(size, size);
    SetWorkingImage(&scaledCharImg);
    cleardevice();
    // 用GetImageHDC获取HDC，StretchBlt实现缩放
    HDC hdcScaled = GetImageHDC(&scaledCharImg);
    HDC hdcTemp = GetImageHDC(&tempCharImg);
    StretchBlt(hdcScaled, 0, 0, size, size,
               hdcTemp, 0, 0, ORIGIN_CHAR_SIZE, ORIGIN_CHAR_SIZE,
               SRCPAINT);
    SetWorkingImage();
    putimage(destX, destY, &scaledCharImg, SRCPAINT);
}

void DrawAsciiString(int startX, int startY, const char* str, int size, IMAGE& g_AsciiFontImg) {
    const int CHAR_SPACING_RATIO = 0;
    int charSpacing = static_cast<int>(size * CHAR_SPACING_RATIO);
    int currentX = startX;
    int currentY = startY;

    for (int i = 0; str[i] != '\0'; i++) {
        DrawSingleAsciiChar(currentX, currentY, str[i], size, g_AsciiFontImg);
        currentX += size + charSpacing;
    }
}

enum GameState {
    SELECT_LEVEL,  // 选关界面
    PLAYING_LEVEL  // 游戏界面
};

int main()
{
    
    const int SLEEP_MS = 10;
    int temp = 0; 
    DWORD lastHitTime = 0; 
    DWORD startTime = 0;
    DWORD startBulletTime = 0;
    const int FADE_DURATION = 3500; 
    bool isDead = false; 
    bool deadMusicPlayed = false;
    bool testmode = false;
    int selectedLevel = 0;
    int frameDelay = 0; 
    const double F_MAX = 3371.0;   //符卡时长 
    const double PI = acos(-1);
    const int WINDOW_WIDTH = 600;
    const int WINDOW_HEIGHT = 700;
    ExMessage msg = { 0 };

	//选关按钮 
    const int BUTTON_WIDTH = 200;
    const int BUTTON_HEIGHT = 60;
    const int BUTTON_SPACING = 40;
    const int BUTTON_X = (WINDOW_WIDTH - BUTTON_WIDTH) / 2;
    const int BUTTON_Y1 = WINDOW_HEIGHT / 2 - BUTTON_HEIGHT - BUTTON_SPACING;
    const int BUTTON_Y2 = WINDOW_HEIGHT / 2;
    const int BUTTON_Y3 = WINDOW_HEIGHT / 2 + BUTTON_HEIGHT + BUTTON_SPACING;

    initgraph(WINDOW_WIDTH + 300, WINDOW_HEIGHT);
    IMAGE background, background1, leftbackground, lives, g_AsciiFontImg;
    loadimage(&background, "assest//eff06b.png", WINDOW_WIDTH, WINDOW_HEIGHT);
    loadimage(&background1, "assest//eff06.png", WINDOW_WIDTH, WINDOW_HEIGHT);
    loadimage(&leftbackground, "assest//leftbackground.png", 300, WINDOW_HEIGHT);
    loadimage(&lives, "assest//lives.png", 30, 30); 
    loadimage(&g_AsciiFontImg, "assest//ascii.png", 256, 96);

    // 灰度背景
    IMAGE grayBackground = background;
    IMAGE grayBackground1 = background1;
    ApplyGrayEffect(&grayBackground);
    ApplyGrayEffect(&grayBackground1);

    // 子弹角度偏移(波与粒的境界)
    const double angleOffsets[] = {8 * PI / 5, 6 * PI /5, 4 * PI / 5, 2 * PI / 5, 0.0};
    const int BULLET_GROUP_COUNT = sizeof(angleOffsets) / sizeof(angleOffsets[0]);
    
    // 子弹参数初始化
    IMAGE Param1Shot,Param2Shot,Param3Shot,Param4Shot;
    BulletParams param1 = {
        &Param1Shot,    // imgShot
        8,          // imgwidth
        16,         // imglength
        8,          // actualsize
        3,        // v
        100.0,      // f
        WINDOW_WIDTH / 2,    // centerX
        WINDOW_HEIGHT / 2,    // centerY
        0.0        // angle
    };//米弹 
    BulletParams param2 = {
        &Param2Shot,    // imgShot
        49,         // imgwidth
        49,         // imglength
        20,         // actualsize
        1.5,        // v
        100.0,      // f
        WINDOW_WIDTH / 2,    // centerX
        WINDOW_HEIGHT / 2,    // centerY
        0.0       // angle
    };//大玉 
    BulletParams param3 = {
        &Param3Shot,    // imgShot
        30,         // imgwidth
        30,         // imglength
        12,         // actualsize
        1,          // v
        100.0,      // f
        WINDOW_WIDTH / 2,    // centerX
        WINDOW_HEIGHT / 2,    // centerY
        0.0       // angle
    };//碟弹 
    BulletParams param4 = {
        &Param4Shot,    // imgShot
        10,         // imgwidth
        10,         // imglength
        5,          // actualsize
        0.5,        // v
        100.0,      // f
        WINDOW_WIDTH / 2,    // centerX
        WINDOW_HEIGHT / 2,    // centerY
        0.0        // angle
    };//小玉 
    
    loadimage(&Param1Shot, "assest//smallshot.png", param1.imgwidth, param1.imglength);
    loadimage(&Param2Shot, "assest//bigshot.png", param2.imgwidth, param2.imglength);
    loadimage(&Param3Shot, "assest//buttershot.png", param3.imgwidth, param3.imglength);
    loadimage(&Param4Shot, "assest//bigshot.png", param4.imgwidth, param4.imglength);

    // 角色初始化
    IMAGE imgChara, imgCharaLeft, imgCharaRight;
    IMAGE imgEnemy, imgEnemyLeft, imgEnemyRight, dead;
    
    Character chara_s = {
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
    
    Character enemy_s ={
        &imgEnemy,           // 默认图像
        &imgEnemyLeft,       // 左移图像
        &imgEnemyRight,      // 右移图像
        144,                 // imgwidth
        84,                  // imglength
        0,                   // actualsize
        3,                   // speed
        0,                   // vx
        0,                   // vy
        300,                 // 初始X坐标
        150,                 // 初始Y坐标
        false,               // showHitPoint
        0                    // direction
    };
    
    loadimage(&imgChara, "assest//player.png", chara_s.imglength, chara_s.imgwidth);
    loadimage(&imgCharaLeft, "assest//playerleft.png", chara_s.imglength, chara_s.imgwidth);
    loadimage(&imgCharaRight, "assest//playerright.png", chara_s.imglength, chara_s.imgwidth);
    loadimage(&imgEnemy, "assest//enemy.png", enemy_s.imglength, enemy_s.imgwidth); 
    loadimage(&imgEnemyLeft, "assest//enemyleft.png", enemy_s.imglength, enemy_s.imgwidth);
    loadimage(&imgEnemyRight, "assest//enemyright.png", enemy_s.imglength, enemy_s.imgwidth);
    loadimage(&dead, "assest//dead.png", enemy_s.imglength, enemy_s.imgwidth);

	const char* stageTexts[3] = {"STAGE 1", "STAGE 2", "STAGE 3"};
	int selectedIndex = 0;
	GameState currentState = SELECT_LEVEL;

    // 游戏状态初始化
    GameStates::ScatterState scatterState;
    GameStates::ButterState butterState;
    GameStates::YuyukoState yuyukoState;
    GameStates::TesseractState tesseractState;

	mciSendStringA("open sounds/select.wav alias select_sound", NULL, 0, NULL);
    mciSendStringA("open sounds/enter.wav alias enter_sound", NULL, 0, NULL);
    mciSendStringA("open sounds/dead.mp3 alias dead_music", NULL, 0, NULL);

    while (true)
    {
        BeginBatchDraw();
        cleardevice();

        // 绘制左侧UI背景
        putimage(600, 0, &leftbackground);
        for(int i=0; i<5-temp; ++i) 
            putimage(700+i*30,123,&lives);
        if(testmode) 
            DrawAsciiString(680, 15, "TEST MODE", 24, g_AsciiFontImg);

        if (currentState == SELECT_LEVEL)
        {
            putimage(0, 0, &grayBackground);
            putimage(0, 0, &grayBackground1, SRCCOPY);
            
            // 绘制标题
            DrawAsciiString(WINDOW_WIDTH/2 - 200, WINDOW_HEIGHT/6, "SELECT YOUR STAGE!", 24, g_AsciiFontImg);
            
            // 处理输入
            while (peekmessage(&msg)) {
                if (msg.message == WM_KEYDOWN) {
                    switch (msg.vkcode) {
                        case VK_UP: {
                            selectedIndex = (selectedIndex - 1 + 3) % 3;
							mciSendStringA("play select_sound from 0", NULL, 0, NULL);
                            break;
                        }
                        case VK_DOWN: {
                            selectedIndex = (selectedIndex + 1) % 3;
							mciSendStringA("play select_sound from 0", NULL, 0, NULL);
                            break;
                        }
                        case 13: // 回车键确认
                        case 'Z': { // Z键确认
	                        mciSendStringA("play enter_sound from 0", NULL, 0, NULL);
                            
                            selectedLevel = selectedIndex + 1;
                            currentState = PLAYING_LEVEL;
                            temp = 0;
                            isDead = false;
                            deadMusicPlayed = false;
                            startTime = GetTickCount();
                            startBulletTime = 0;
                            param1.f = 0.0;
                            lastHitTime = 0;
							chara_s.vx = 0;
							chara_s.vy = 0;
							chara_s.direction = 0; 
                            
                            mciSendString("stop bgm", NULL, 0, NULL);
							mciSendString("close bgm", NULL, 0, NULL);
							mciSendString("stop dead_music", NULL, 0, NULL);
							mciSendString("close dead_music", NULL, 0, NULL);
                 	    	mciSendStringA("open sounds/background_music.mp3 alias bgm", NULL, 0, NULL);
            	            mciSendStringA("play bgm repeat", NULL, 0, NULL);
                            
                            // 重置角色位置
                            enemy_s.cX = 300;
                            enemy_s.cY = 150;
                            chara_s.cX = WINDOW_WIDTH / 2;
                            chara_s.cY = WINDOW_HEIGHT - 32;
                            
                            // 重置子弹状态
                            GameStates::ResetScatterState(scatterState);
                            GameStates::ResetButterState(butterState);
                            GameStates::ResetYuyukoState(yuyukoState);
                            GameStates::ResetTesseractState(tesseractState);
                            break;
                        }
                    }
                }
            }
            
            // 绘制关卡选择
            const int BASE_SIZE = 24;
            const int SELECTED_SIZE = 32;
            const int SPACING = 60;
            int startY = WINDOW_HEIGHT / 3;
            
            for (int i = 0; i < 3; i++) {
                int textSize = (i == selectedIndex) ? SELECTED_SIZE : BASE_SIZE;
                int textWidth = strlen(stageTexts[i]) * textSize;
                int textX = WINDOW_WIDTH / 2 - textWidth / 2;
                int textY = startY + i * SPACING;
                DrawAsciiString(textX, textY, stageTexts[i], textSize, g_AsciiFontImg);
            }
        }
        else if (currentState == PLAYING_LEVEL)
        {
	        putimage(0, 0, &background, SRCCOPY);
	        putimage(0, 0, &background1, SRCPAINT);
            
            DWORD currentTime = GetTickCount();
            if (currentTime - startTime >= 1000 && !isDead) {
                if (startBulletTime == 0) 
                    startBulletTime = currentTime;

                double progress = static_cast<double>(currentTime - startBulletTime) / FADE_DURATION;
                progress = progress > 1.0 ? 1.0 : progress;
                progress = progress < 0.0 ? 0.0 : progress;

                // 根据关卡绘制子弹
                switch(selectedLevel){
                    case 1:
                        for (int i = 0; i < BULLET_GROUP_COUNT; ++i) {
                            param1.angleOffset = angleOffsets[i];
                            DrawBullet_yokary(param1, chara_s, &temp, &lastHitTime, progress, isDead, i);
                        }      
                        DrawBullet_butter(butterState, param3, chara_s, enemy_s, &temp, &lastHitTime, isDead);
                        break;
                    case 2:
                        DrawDemo(tesseractState, param4, chara_s, &temp, &lastHitTime, isDead);
                        DrawBullet_butter(butterState, param3, chara_s, enemy_s, &temp, &lastHitTime, isDead);
                        DrawBullet_yuyuko(yuyukoState, param2, chara_s, &temp, &lastHitTime, progress, isDead);
                        break;
                    case 3:
                        DrawScatterBullet(scatterState, param4, chara_s, &temp, &lastHitTime, progress, isDead);
                        DrawBullet_butter(butterState, param3, chara_s, enemy_s, &temp, &lastHitTime, isDead);
                        break;
                }
            }
            
            // 绘制敌人
            if (!isDead){
                IMAGE* enemyCurrentImg = enemy_s.imgChara; 
                if (enemy_s.direction == 1) {
                    enemyCurrentImg = enemy_s.imgCharaLeft;
                } else if (enemy_s.direction == 2) {
                    enemyCurrentImg = enemy_s.imgCharaRight;
                }
                putimage(enemy_s.cX - enemy_s.imglength/2, enemy_s.cY - enemy_s.imgwidth/2, enemyCurrentImg, SRCPAINT);
            }
            
            // 绘制玩家
            if (!isDead && (currentTime - lastHitTime >= 1514 || (currentTime % 200) < 100)) {
                IMAGE* currentImg = chara_s.imgChara;
                if (chara_s.direction == 1) { 
                    currentImg = chara_s.imgCharaLeft;
                } else if (chara_s.direction == 2) { 
                    currentImg = chara_s.imgCharaRight;
                }
                putimage(chara_s.cX - chara_s.imglength/2, chara_s.cY - chara_s.imgwidth/2, currentImg, SRCPAINT);
                if (chara_s.showHitPoint) {
                    setfillcolor(WHITE);
                    setlinecolor(RED);
                    fillcircle(chara_s.cX, chara_s.cY, chara_s.actualsize);
                }
            }
            
            // 死亡判断
            if (temp >= 5) {
                isDead = true;
                ApplyGrayEffect(&grayBackground);
                ApplyGrayEffect(&grayBackground1);
                if (!deadMusicPlayed) {
                	mciSendString("stop bgm", NULL, 0, NULL);
					mciSendStringA("play dead_music from 0", NULL, 0, NULL); 
                    deadMusicPlayed = true;
                }
            }
            else if (param1.f > F_MAX){
            	isDead = true;
            	mciSendString("close bgm", NULL, 0, NULL);
            	currentState = SELECT_LEVEL;
			} 
            
            if (isDead) {
                settextcolor(WHITE);
                setbkmode(TRANSPARENT);
                DrawAsciiString(WINDOW_WIDTH/2 - 105, WINDOW_HEIGHT/2, "Game Over", 24, g_AsciiFontImg);
                putimage(WINDOW_WIDTH/2 - enemy_s.imgwidth/2, 100, &dead, SRCPAINT);
                while (peekmessage(&msg)) {
                    if (msg.message == WM_KEYDOWN && msg.vkcode == 'X'){
                        currentState = SELECT_LEVEL;
                        mciSendString("stop deadbgm", NULL, 0, NULL);
                        mciSendString("close deadbgm", NULL, 0, NULL);
                        deadMusicPlayed = false;
                    }
                }
            }
            
            // 玩家控制
            while (peekmessage(&msg, EX_KEY)) {
                if (msg.message == WM_KEYDOWN) {
                    switch (msg.vkcode) {
                        case VK_UP:    chara_s.vy = -1.5; break;
                        case VK_DOWN:  chara_s.vy = 2;  break;
                        case VK_LEFT:  
                            chara_s.vx = -1.5; 
                            chara_s.direction = 1;
                            break;
                        case VK_RIGHT: 
                            chara_s.vx = 1.75;  
                            chara_s.direction = 2; 
                            break;
                        case VK_SHIFT: 
                            chara_s.speed = 1; 
                            chara_s.showHitPoint = true;
                            break;
                        case 'T': 
                            testmode = !testmode;
							mciSendStringA("play select_sound from 0", NULL, 0, NULL); 
                            break;
                        case 'X': 
                            currentState = SELECT_LEVEL;
                            mciSendStringA("play enter_sound from 0", NULL, 0, NULL);
                            mciSendString("stop bgm", NULL, 0, NULL);
                            mciSendString("close bgm", NULL, 0, NULL);
                            break;
                    }
                } else if (msg.message == WM_KEYUP) {
                    switch (msg.vkcode) {
                        case VK_UP:    chara_s.vy = 0; break;
                        case VK_DOWN:  chara_s.vy = 0; break;
                        case VK_LEFT:  
                            chara_s.vx = 0; 
                            if (chara_s.direction == 1) chara_s.direction = 0;
                            break;
                        case VK_RIGHT: 
                            chara_s.vx = 0; 
                            if (chara_s.direction == 2) chara_s.direction = 0;
                            break;
                        case VK_SHIFT: 
                            chara_s.speed = 2; 
                            chara_s.showHitPoint = false;
                            break;
                    }
                }
            }
    
            // 角色移动和边界检测
            chara_s.cX += chara_s.speed * chara_s.vx;
            chara_s.cY += chara_s.speed * chara_s.vy;  
            if (chara_s.cX - chara_s.imglength/2 < 0) 
                chara_s.cX = chara_s.imglength/2;
            if (chara_s.cX + chara_s.imglength/2 > WINDOW_WIDTH) 
                chara_s.cX = WINDOW_WIDTH - chara_s.imglength/2;
            if (chara_s.cY - chara_s.imgwidth/2 < 0) 
                chara_s.cY = chara_s.imgwidth/2;
            if (chara_s.cY + chara_s.imgwidth/2 > WINDOW_HEIGHT) 
                chara_s.cY = WINDOW_HEIGHT - chara_s.imgwidth/2;
            
            enemy_s.actualsize = enemy_s.cX;
            if (!isDead) UpdateEnemyMovement(enemy_s);
            
            // 作弊模式
            if (testmode) temp = 0;
                
            // 显示进度
            char tempstr[20];
            if (!isDead) {
                sprintf(tempstr, "%d/%d", (int)param1.f, (int)F_MAX);
            }
            DrawAsciiString(WINDOW_WIDTH+100, WINDOW_HEIGHT/2-95, tempstr, 21, g_AsciiFontImg);
        }

        EndBatchDraw();
        Sleep(SLEEP_MS);
        param1.f += 1.0;
        
        if (GetAsyncKeyState(VK_ESCAPE) & 0x8000) {
            break;
        }
    }

    mciSendString("close bgm", NULL, 0, NULL);
    closegraph();
    return 0;
}
