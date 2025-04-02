#include <iostream>
#include <vector>
#include <windows.h>
#include <conio.h>
#include <time.h>
#include <fstream>
#include <string>
#include <queue>

using namespace std;
const int cNum = 48;

HANDLE hConsole;
HANDLE hStdOut = GetStdHandle(STD_OUTPUT_HANDLE);
COORD position;

void GotoXY(int X, int Y) {
    COORD pos = { X, Y };
    SetConsoleCursorPosition(hStdOut, pos);
}

enum class TypeCell {
    dead,
    alive,
    predator,
    barrier,
    portal
};

struct Cell {
    TypeCell typeCell = TypeCell::dead;
    static inline char symbols[5] = { '.', '#', '*', 'N', '@'};

    friend ostream& operator<<(ostream& out, const Cell& cell) {
        out << Cell::symbols[(int)cell.typeCell];
        return out;
    }
};

struct Predator {
    pair<int, int> pos;
    int livetime = 0;
    Predator(int i, int j) {
        pos.first = i;
        pos.second = j;
    }
};

struct Portal {
    pair<int, int> pos;
    pair<int, int> friend_pos;
    int number;
    Portal(int i, int j) {
        pos.first = i;
        pos.second = j;
        //number = n;
    }
};

int PortalCounter = 0;

class Field
{
private:
    int n; // число строк
    int m; // число столбцов
public:
    vector<vector<Cell>> cells;
    vector<Predator> Predators;
    vector<Portal> Portals;
    Field() {};
    Field(int n, int m, float barrierProportion = 0.0f, float aliveProportion = 0.3f, float predatorProportion = 0.1f, float portalProportion = 0.0f) : n(n), m(m) {
        srand(0);
        // случайное поле n на m, с шансом появления живой клетки = aliveProportion
        cells = vector<vector<Cell>>(n);
        for (int i = 0; i < n; i++) {
            cells[i] = vector<Cell>(m);
            for (int j = 0; j < m; j++) {
                Cell cell;
                float prop = (rand() % 10000) / 10000.f;
                if (prop <= portalProportion) {
                    if (j == m - 1)
                        cell.typeCell = TypeCell::dead;
                    else {
                        cell.typeCell = TypeCell::portal;
                        Portals.push_back(Portal(i, j));
                        if (Portals.size() % 2 == 0) {
                            Portals[Portals.size() - 2].friend_pos.first = Portals[Portals.size() - 1].pos.first;
                            Portals[Portals.size() - 2].friend_pos.second = Portals[Portals.size() - 1].pos.second;
                            Portals[Portals.size() - 1].friend_pos.first = Portals[Portals.size() - 2].pos.first;
                            Portals[Portals.size() - 1].friend_pos.second = Portals[Portals.size() - 2].pos.second;
                            Portals[Portals.size() - 2].number = PortalCounter++;
                            Portals[Portals.size() - 1].number = Portals[Portals.size() - 2].number;
                        }
                        cells[i][j] = cell;
                        j++;
                    }
                }
                else if (prop <= predatorProportion) {
                    cell.typeCell = TypeCell::predator;
                    Predators.push_back(Predator(i, j));
                    cells[i][j] = cell;
                }
                else if (prop <= aliveProportion) {
                    cell.typeCell = TypeCell::alive;
                    cells[i][j] = cell;
                }
                else if (prop <= barrierProportion) {
                    cell.typeCell = TypeCell::barrier;
                    cells[i][j] = cell;
                }
                else {
                    cell.typeCell = TypeCell::dead;
                    cells[i][j] = cell;
                }
            }
        }
        if (Portals.size() % 2 == 1) {
            int i = Portals[Portals.size() - 1].pos.first;
            int j = Portals[Portals.size() - 1].pos.second;
            Portals.pop_back();
            cells[i][j].typeCell = TypeCell::dead;
        }
    }

    int sizeh() const {
        return this->n;
    }

    int sizew() const {
        return this->m;
    }

    friend ostream& operator<<(ostream& out, const Field& field) {
        out << endl << field.n << ' ' << field.m << endl;
        for (int i = 0; i < field.n; i++) {
            for (int j = 0; j < field.m; j++) {
                out << field.cells[i][j];
                if (field.cells[i][j].typeCell == TypeCell::portal)
                    for (int x = 0; x < field.Portals.size(); x++)
                        if (field.Portals[x].pos.first == i && field.Portals[x].pos.second == j) {
                            cout << field.Portals[x].number;
                            j++;
                        }
            }
            out << endl;
        }
        return out;
    }

    vector<Cell>& operator[](int i) {
        return cells[i];
    }

    Field& operator=(const Field& f) {
        n = f.n;
        m = f.m;
        cells.clear();
        cells = vector<vector<Cell>>(n);
        for (int i = 0; i < n; i++) {
            cells[i] = vector<Cell>(m);
            for (int j = 0; j < m; j++) {
                cells[i][j] = f.cells[i][j];
            }
        }
        Predators.clear();
        for (int i = 0; i < f.Predators.size(); i++)
            Predators.push_back(f.Predators[i]);
        Portals.clear();
        for (int i = 0; i < f.Portals.size(); i++)
            Portals.push_back(f.Portals[i]);
        return *this;
    }
    // возвращает сколько клеток типа type находится в радиусе=radius, вокруг клетки posX, posY
    int getNum(int posX, int posY, TypeCell type = TypeCell::alive, int radius = 1) const {
        int res = 0, i, j, I, J;
        i = posY - radius;
        j = posX - radius;
        I = posY + radius;
        J = posX + radius;
        while (i <= I) {
            if (i >= 0 && i <= n - 1) {
                while (j <= J) {
                    if (j >= 0 && j <= m - 1) {
                        if ((posY == i && posX == j)) {
                            j++;
                            continue;
                        }
                        if (cells[i][j].typeCell != type) {
                            j++;
                        }
                        else {
                            j++;
                            res++;
                        }
                    }
                    else
                        j++;
                }
            }
            i++;
            j = posX - radius;
        }
        return res;
    }

    bool chekPortal(int i, int j) {
        if (cells[i][j].typeCell == TypeCell::portal)
            return false;
        else return true;
    }

    bool checkPos(int i, int j) {
        if (i >= sizeh() || i < 0 || j >= sizew() || j < 0 || (j > 0 && cells[i][j - 1].typeCell == TypeCell::portal) || cells[i][j].typeCell == TypeCell::barrier)
            return false;
        else return true;
    }

    pair<int, int> getMove(int i, int j, bool ifportal = false) {
        pair<int, int> offsets[8] = { {-1, -1}, {0, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0} };
        if (ifportal)
            pair<int, int> offsets[7] = { {-1, -1}, {1, -1}, {1, 0}, {1, 1}, {0, 1}, {-1, 1}, {-1, 0} };
        int n = sizeof(offsets) / sizeof(offsets[0]);
        int randStep = time(0) % n;
        int count = 0;
        pair<int, int> nextPos = { i + offsets[randStep].first, j + offsets[randStep].second };
        if (!ifportal) {
            while (!checkPos(nextPos.first, nextPos.second) && count < n) {
                randStep = (randStep + 1) % n;
                count++;
                nextPos = { i + offsets[randStep].first, j + offsets[randStep].second };
            }
        }
        else
            while ((!checkPos(nextPos.first, nextPos.second) && count < n) || !chekPortal(nextPos.first, nextPos.second)) {
                randStep = (randStep + 1) % n;
                count++;
                nextPos = { i + offsets[randStep].first, j + offsets[randStep].second };
            }
        if (count == n)
            return { i, j };
        return nextPos;
    }
};

void updateAllConsole(Field maingame) {
    cout << maingame;
    SetConsoleCursorPosition(hConsole, position);
}

struct iGame {
    int n = 0;
    int m = 0;

    int seed = 0; // случайная величина для генератора
    double aliveProportion = 0.0;  // вероятность того, что клетка живая
    double predatorProportion = 0.0;  // вероятность того, что клетка хищная
    double barrierProportion = 0.4;  // вероятность того, что клетка хищная
    double predatorDivPropotion = 0.3; // вероятность деления хищной клетки
    double ldProportion = 0.01; // вероятность возникновения живой клетки на месте мертвой
    int dimension = 2; // размерность

    int radius = 1; // радиус проверки, граница включена
    int loneliness = 1; // с этого числа и меньше клетки умирают от одиночества
    int birth_start = 3; // с этого числа и до birth_end появляется живая клетка
    int birth_end = 2;
    int overpopulation = 4; // с этого числа и дальше клетки погибают от перенаселения
    int predator_jump = 1; // число, на которое смещается хищная клетка за один ход
    int predator_livetime = 20;

    virtual void runGame(int numIt) = 0;
    virtual bool save(const string& filename) const = 0;
    virtual bool load(const string& filename) = 0;
};

Field Predator_move(int i, int j, Field field, int x) {
    pair<int, int> newPos;
    do newPos = field.getMove(i, j);
    while (newPos.second > 0 && field.cells[newPos.first][newPos.second - 1].typeCell == TypeCell::portal);
    int i1 = newPos.first, j1 = newPos.second;//i1, j1 - новые координаты хищника
    if (j1 > 0 && field[i1][j1 - 1].typeCell == TypeCell::portal)
        j1 = j1 - 1;
    if (field[i1][j1].typeCell == TypeCell::alive) {
        field[i1][j1].typeCell = TypeCell::predator;
        field.Predators[x].pos.first = i1;
        field.Predators[x].pos.second = j1;
        field.Predators[x].livetime = 0;
        double prop = (rand() % 10000) / 10000.f;
        if (prop <= 0.3) {
            field[i][j].typeCell = TypeCell::predator;//если хищник поделился
            field.Predators.push_back(Predator(i, j));
        }
        else
            field[i][j].typeCell = TypeCell::dead;
    }
    else if (field[i1][j1].typeCell == TypeCell::dead) {
        field[i][j].typeCell = TypeCell::dead;
        if (field.Predators[x].livetime == 20) {
            swap(field.Predators[x], field.Predators[field.Predators.size() - 1]);
            field.Predators.pop_back();
        }
        else {
            field[i1][j1].typeCell = TypeCell::predator;
            field.Predators[x].livetime++;
            field.Predators[x].pos.first = i1;
            field.Predators[x].pos.second = j1;
        }
    }
    else if (field[i1][j1].typeCell == TypeCell::portal)
        for (int y = 0; y < field.Portals.size(); y++) {
            if (field.Portals[y].pos.first == i1 && field.Portals[y].pos.second == j1) {
                pair<int, int> newpos = field.getMove(field.Portals[y].friend_pos.first, field.Portals[y].friend_pos.second, true);
                if (field[newpos.first][newpos.second].typeCell == TypeCell::alive) {
                    field.Predators[x].livetime = 0;
                    field[newpos.first][newpos.second].typeCell = TypeCell::predator;
                    field.Predators[x].pos.first = newpos.first;
                    field.Predators[x].pos.second = newpos.second;
                }
                else if (field[newpos.first][newpos.second].typeCell == TypeCell::dead) {
                    if (field.Predators[x].livetime == 20) {
                        field[newpos.first][newpos.second].typeCell = TypeCell::dead;
                        swap(field.Predators[x], field.Predators[field.Predators.size() - 1]);
                        field.Predators.pop_back();
                    }
                    else {
                        field.Predators[x].livetime++;
                        field[newpos.first][newpos.second].typeCell = TypeCell::predator;
                        field.Predators[x].pos.first = newpos.first;
                        field.Predators[x].pos.second = newpos.second;
                    }
                }
                field[i][j].typeCell = TypeCell::dead;
                break;
            }
        }
    else if (i1 == i && j1 == j) {
        if (field.Predators[x].livetime == 20) {
            swap(field.Predators[x], field.Predators[field.Predators.size() - 1]);
            field.Predators.pop_back();
            field[i][j].typeCell = TypeCell::dead;
        }
        else field.Predators[x].livetime++;
    }
    else { //если хищник попал на хищника
        field[i][j].typeCell = TypeCell::dead;
        swap(field.Predators[x], field.Predators[field.Predators.size() - 1]);
        field.Predators.pop_back();
        for (int y = 0; y < field.Predators.size(); y++) {
            if (field.Predators[y].pos.first == i1 && field.Predators[y].pos.second == j1)
                field.Predators[y].livetime = 0;
            break;
        }
    }
    return field;
}

struct Game : public iGame {
    Field maingame;

    Game() {}

    Game(int n, int m) {
        Field a(n, m);
        maingame = a;
        cout << maingame;
    }

    void runstep()
    {
        Field a = maingame;
        Field field(a.sizeh(), a.sizew());
        for (int i = 0; i < a.sizeh(); i++) {
            for (int j = 0; j < a.sizew(); j++) {
                if (a[i][j].typeCell == TypeCell::alive) {
                    if (a.getNum(j, i) > iGame::loneliness && a.getNum(j, i) < iGame::overpopulation)
                        field[i][j].typeCell = TypeCell::alive;
                    else
                        field[i][j].typeCell = TypeCell::dead;
                }
                else if (a[i][j].typeCell == TypeCell::dead) {
                    if (j>0 && a[i][j - 1].typeCell == TypeCell::portal)
                        a[i][j].typeCell == TypeCell::dead;
                    else if (a.getNum(j, i) <= iGame::birth_start && a.getNum(j, i) >= iGame::birth_end)
                        field[i][j].typeCell = TypeCell::alive;
                    else {
                        float prop = (rand() % 10000) / 10000.f;
                        if (prop <= ldProportion)
                            field[i][j].typeCell = TypeCell::alive;
                        else
                            field[i][j].typeCell = TypeCell::dead;
                    }
                }
                else if (a[i][j].typeCell == TypeCell::barrier)
                    field[i][j].typeCell = TypeCell::barrier;
                else if (a[i][j].typeCell == TypeCell::portal)
                    field[i][j].typeCell = TypeCell::portal;
            }
        }
        field.Predators = a.Predators;
        field.Portals = a.Portals;
        for (int x = 0; x < field.Predators.size(); x++) {
            int i = field.Predators[x].pos.first, j = field.Predators[x].pos.second;
            field = Predator_move(i, j, field, x);
        }
        maingame = field;
    }

    bool save(const string& results) const override {
        ofstream out;
        out.open("results.txt");
        if (out.is_open()) {
            out << maingame.sizeh() << ' ' << maingame.sizew() << ' ' << "'#' - alive" << ' ' << "'.' - dead" << "'*' - predator" << "'&' - barrier" << "'@' - portal" << endl;
            for (int i = 0; i < maingame.sizeh(); i++) {
                for (int j = 0; j < maingame.sizew(); j++)
                    out << maingame.cells[i][j];
                out << endl;
            }
        }
        out.close();
        return true;
    }

    bool load(const string& results) override {
        ifstream in("results.txt");
        if (in.is_open()) {
            string first, strnum;
            getline(in, first);
            int rows = 0, col = 0, i = 0, number = 0;
            for (i; first[i] != ' '; i++)
                strnum.push_back(first[i]);
            number = stoi(strnum);
            rows = number;
            number = 0;
            i++;
            strnum = "";
            for (i; first[i] != ' '; i++)
                strnum.push_back(first[i]);
            number = stoi(strnum);
            col = number;
            char alive = '#';
            char dead = '.';
            char predator = '*';
            char barrier = '&';
            char portal = '@';
            Field temp(rows, col, 0, 0, 0);
            for (int i = 0; i < rows; i++) {
                getline(in, first);
                for (int j = 0; j < col; j++) {
                    if (first[j] == alive)
                        temp[i][j].typeCell = TypeCell::alive;
                    else if (first[j] == dead)
                        temp[i][j].typeCell = TypeCell::dead;
                    else if (first[j] == predator) {
                        temp[i][j].typeCell = TypeCell::predator;
                        temp.Predators.push_back(Predator(i, j));
                    }
                    else if (first[j] == portal) {
                        temp[i][j].typeCell = TypeCell::portal;
                        temp.Portals.push_back(Portal(i, j));
                        temp.Portals[temp.Portals.size() - 1].number = first[j + 1] - cNum;
                        j++;
                        if (temp.Portals.size() % 2 == 0) {
                            temp.Portals[temp.Portals.size() - 2].friend_pos.first = temp.Portals[temp.Portals.size() - 1].pos.first;
                            temp.Portals[temp.Portals.size() - 2].friend_pos.second = temp.Portals[temp.Portals.size() - 1].pos.second;
                            temp.Portals[temp.Portals.size() - 1].friend_pos.first = temp.Portals[temp.Portals.size() - 2].pos.first;
                            temp.Portals[temp.Portals.size() - 1].friend_pos.second = temp.Portals[temp.Portals.size() - 2].pos.second;
                        }
                    }
                    else
                        temp[i][j].typeCell = TypeCell::barrier;
                }
            }
            maingame = temp;
            cout << endl << "Loaded field:" << endl << maingame << endl << endl;
        }
        return false;
    }

    void runGame(int numIt) override {
        cout << maingame;
        for (int i = 0; i < numIt; i++) {
            system("cls");
            if (_kbhit()) {
                int button = _getch();
                if (button) {
                    if (button == 27)//esc
                        exit(0);
                    else if ((char)button == 's')
                        save("results.txt");
                    else if ((char)button == 'l') {
                        load("results.txt");
                    }
                }
            }
            this->runstep();
            cout << maingame << endl << endl;
            cout << i;
            Sleep(1000);
        }
        //save("results.txt");
    }
    friend Field;
};

struct newGame : public Game {
    newGame(int n, int m) {
        Field a(n, m, 0.4, 0.3, 0.1, 0.05);
        maingame = a;
        cout << maingame;
    }
};


pair<vector<vector<int>>, vector<vector<int>>> lists(Field field) {
    vector<vector<int>> temp(field.sizeh());
    int count = 0;
    int x = 0, y = 0;
    for (int i = 0; i < field.sizeh(); i++) {
        temp[i] = vector<int>(field.sizew());
        for (int j = 0; j < field.sizew(); j++) {
            if (field[i][j].typeCell != TypeCell::barrier) {
                temp[i][j] = count;
                count++;
            }
            else
                temp[i][j] = -1;
        }
    }
    vector<vector<int>> res;// (count);
    int h = 0;
    for (int i = 0; i < field.sizeh(); i++)
        for (int j = 0; j < field.sizew(); j++)
            if (temp[i][j] != -1)
                for (int k = -1; k < 2; k++)
                    for (int l = -1; l < 2; l++)
                        if (!(k == 0 && l == 0))
                            if (field.checkPos(i + k, j + l))
                                    res[temp[i][j]].push_back(temp[i + k][j + l]);
    pair<vector<vector<int>>, vector<vector<int>>> t(temp, res);
    return t;
}

pair<int, int> numPos(vector<vector<int>> a, int number, Field field) {
    for (int i = 0; i < field.sizeh(); i++)
        for (int j = 0; j < field.sizew(); j++)
            if (a[i][j] == number)
                return make_pair(i, j);
}

pair<int, int> bfs(vector<vector<int>>& graph, int start, Field field, vector<vector<int>> a) {
    vector<bool> visited(graph.size());
    vector<vector<int>> dist(graph.size());
    queue<int> q;
    int x0 = numPos(a, start, field).second;
    int y0 = numPos(a, start, field).first;
    visited[start] = true;
    dist[start].push_back(start);
    q.push(start);
    while (!q.empty()) {
        int tek = q.front();
        q.pop();
        for (int t : graph[tek]) {
            int y = numPos(a, t, field).first;
            int x = numPos(a, t, field).second;
            if (!visited[t]) {
                visited[t] = true;
                for (int y : dist[tek])
                    dist[t].push_back(y);
                dist[t].push_back(t);
                q.push(t);
            }
            if (field[y][x].typeCell == TypeCell::alive) {
                int r = dist[t][1];
                int y1 = numPos(a, r, field).first;
                int x1 = numPos(a, r, field).second;
                pair<int, int> res;
                res.first = y1;
                res.second = x1;
                return res;
            }
        }
    }
    return make_pair(y0, x0);
}

int x = 0, y = 0;

struct newPredator : protected Game{
private:
    vector<vector<int>> graph;
    vector<vector<int>> smej;
    Field field;
    void runstep(Field a) {
        for (int i = 0; i < field.sizeh(); i++)
            for (int j = 0; j < field.sizew(); j++)
                if (field[i][j].typeCell == TypeCell::predator)
                    for(int x=0; x<field.Predators.size(); x++)
                        if (field.Predators[x].pos.first == i && field.Predators[x].pos.second == j) {
                            /*if (field.Predators[x].livetime >= predator_livetime) {
                                a[i][j].typeCell = TypeCell::dead;
                            }
                            else {*/
                                pair<int, int> move = bfs(smej, graph[i][j], a, graph);

                                //int h = unit.lifetime;
                                Predator_move(i, j, a, x);
                            //}
                        }
        field = a;
        updateAllConsole(field);
        for (int i = 0; i < field.sizeh(); i++)
            for (int j = 0; j < field.sizew(); j++) {
                if (a[i][j].typeCell == TypeCell::dead) {
                    if (a.getNum(j, i) >= birth_start && a.getNum(j, i) <= birth_end)
                        field[i][j].typeCell = TypeCell::alive;
                    else if (rand() % 100 <= ldProportion)
                        field[i][j].typeCell = TypeCell::alive;
                    else
                        field[i][j].typeCell = TypeCell::dead;
                }
                else if (field[i][j].typeCell == TypeCell::alive) {
                    if (a.getNum(j, i) > loneliness && a.getNum(j, i) < overpopulation)
                        field[i][j].typeCell = TypeCell::alive;
                    else
                        field[i][j].typeCell = TypeCell::dead;
                }
            }
    }
public:
    newPredator(int n, int m) {
        field = Field(n, m, 0.4, 0.3, 0.1, 0.0);
        /*for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                if (rand() % 10 > barrierProportion * 10)
                    field[i][j].typeCell = TypeCell::barrier;*/
        updateAllConsole(field);
        //Sleep(1500);
        graph = lists(field).first;
        smej = lists(field).second;
    };
    void runGame(int numIt) override {
        int i = 0;
        bool flag = 0, spaceflag = 0, mainflag = 0;
        while (i < numIt){
            if (_kbhit()){
                char button = _getch();
                if (button == 32 && spaceflag == 0){
                    while (spaceflag == 0){
                        button = 0;
                        while (button != 32){
                            if (_kbhit())
                                button = _getch();
                            if (button == 'd') x = (x + 1 + field.sizew()) % field.sizew();
                            else if (button == 'a') x = (x - 1 + field.sizew()) % field.sizew();
                            else if (button == 'w') y = (y - 1 + field.sizeh()) % field.sizeh();
                            else if (button == 'x') y = (y + 1 + field.sizeh()) % field.sizeh();
                            else if (button == 'r'){
                                if (field[y][x].typeCell == TypeCell::dead){
                                    field[y][x].typeCell = TypeCell::alive;
                                    putchar('#');
                                }
                                else if (field[y][x].typeCell == TypeCell::alive){
                                    field[y][x].typeCell = TypeCell::predator;
                                    field.Predators.push_back(Predator(x, y));
                                    putchar('*');
                                }
                                else if (field[y][x].typeCell == TypeCell::predator){
                                    field[y][x].typeCell = TypeCell::barrier;
                                    putchar('N');
                                }
                                else if (field[y][x].typeCell == TypeCell::barrier){
                                    field[y][x].typeCell = TypeCell::dead;
                                    putchar('.');
                                }
                            }
                            position.X = x;
                            position.Y = y;
                            SetConsoleCursorPosition(hConsole, position);
                            if (button == 27){
                                position.X = 0;
                                position.Y = field.sizeh() + 1;
                                SetConsoleCursorPosition(hConsole, position);
                                cout << "wasted";
                                Sleep(1000);
                                mainflag = 1;
                                spaceflag = 1;
                                break;
                            }
                            else if ((char)button == 's'){
                                save("results.txt");
                                position.X = 0;
                                position.Y = field.sizeh() + 1;
                                SetConsoleCursorPosition(hConsole, position);
                                cout << "saved";
                                Sleep(1000);
                                position.X = 0;
                                position.Y = field.sizeh() + 1;
                                SetConsoleCursorPosition(hConsole, position);
                                cout << "      ";
                                position.X = 0;
                                position.Y = 0;
                                SetConsoleCursorPosition(hConsole, position);
                                button = 0;
                            }
                            else if ((char)button == 'l'){
                                load("saves.txt");
                                position.X = 0;
                                position.Y = 0;
                                SetConsoleCursorPosition(hConsole, position);
                                updateAllConsole(field);
                                button = 0;
                            }
                            if (button == 32){
                                spaceflag = 1;
                                position.X = 0;
                                position.Y = 0;
                                SetConsoleCursorPosition(hConsole, position);
                                continue;
                            }
                            button = 0;
                        }
                    }
                }
            }
            else{
                runstep(field);
                i++;
                updateAllConsole(field);
                Sleep(1500);
                spaceflag = 0;
            }
            if (mainflag == 1)
                break;
        }
    }
};

class Zmey {
public:
    COORD* t;
    int FCount;
};

enum{  //направление змейки
    Left,
    Up,
    Right,
    Down
};

class ZGame {
public:
    Zmey zmey;
    COORD food;
    int dx, dy;
    int pause;
    int dir;

    enum{
        end,
        kray,
        plus,
        move
    };

    void PlusFood() { // Распределение еды
        int i = 0, x, y;
        int n = zmey.FCount;
        while (i < n) {
            x = rand() % 56 + 3; //координаты еды
            y = rand() % 19 + 3;
            for (i = 0; i < n; i++)
                if (x == zmey.t[i].X && y == zmey.t[i].Y)
                    break;
        }
        food.X = x;
        food.Y = y;
        SetConsoleCursorPosition(hConsole, food);
        printf("%c", '#');
    }

    void sGame() {
        system("cls");
        zmey.FCount = 3;
        zmey.t = new COORD[3];
        for (int i = 0; i < 3; i++) {
            zmey.t[i].X = 20 + i;
            zmey.t[i].Y = 20;
        }
        dx = 1;
        dy = 0;
        pause = 100;
        PlusFood();
    }

    void Kray() {
        GotoXY(0, 0);
        for (int i = 0; i < 62; i++)
            printf("*");
        GotoXY(0, 24);
        for (int i = 0; i < 62; i++)
            printf("*");
        GotoXY(0, 1);
        for (int i = 0; i < 24; i++)
            cout << "*" << endl;
        for (int i = 1; i < 24; i++) {
            GotoXY(61, i);
            cout << "*" << endl;
        }
    }

    int Move() {
        int& n = zmey.FCount;
        COORD head = zmey.t[n - 1];
        COORD tail = zmey.t[0];
        COORD next;
        next.X = head.X + dx;
        next.Y = head.Y + dy;
        if (next.X < 1 || next.Y < 1 || next.X > 60 || next.Y > 23)
            return kray;
        if (n > 4)
            for (int i = 0; i < n; i++)
                if (next.X == zmey.t[i].X && next.Y == zmey.t[i].Y)
                    return end;
        if (next.X == food.X && next.Y == food.Y) {
            COORD* temp = new COORD[++n];
            for (int i = 0; i < n; i++)
                temp[i] = zmey.t[i];
            temp[n - 1] = next;
            delete[] zmey.t;
            zmey.t = temp;
            SetConsoleCursorPosition(hConsole, head);
            cout << "*";
            SetConsoleCursorPosition(hConsole, next);
            cout << "*";
            PlusFood();
            return plus;
        }
        for (int i = 0; i < n - 1; i++)
            zmey.t[i] = zmey.t[i + 1];
        zmey.t[n - 1] = next;
        SetConsoleCursorPosition(hConsole, tail);
        cout << " ";
        SetConsoleCursorPosition(hConsole, head);
        cout << "*";
        SetConsoleCursorPosition(hConsole, next);
        cout << "@";
        return move;
    }
};

int main() {
    //cout << "Press 'esc' to exit" << endl << "\t's' - to save" << endl << "\t'l' - to load" << endl << "\t'enter' - to start game";
    ///*Game test(5, 30);
    //if (_getch() == 13)
    //    test.runGame(23);*/
    ///*newGame test(10, 20);
    //if (_getch() == 13)
    //    test.runGame(30);*/

    //newPredator test(3, 5);
    //if (_getch() == 13)
    //    test.runGame(30);
    cout << "Use arrows on the keyboard to control snake." << endl << "Press 'enter' - to start game";
    if (_getch() == 13) {
        hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
        int key = 0, count = 0;
        bool Pause = false;
        ZGame game;
        game.sGame();
        game.Kray();
        srand(time(0));
        bool pause = false;
        while (key != 27) {
            while (!_kbhit()) {
                if (Pause == true) {
                    Sleep(1);
                    continue;
                }
                switch (game.Move()) {
                case game.plus:
                    ++count;
                    game.pause -= 1;
                    if (count == 50) {
                        system("cls");
                        cout << "Victory!";
                        return(0);
                    }
                    break;
                case game.kray:
                case game.end:
                    system("cls");
                    cout << "Game over.";
                    return(0);
                }
                Sleep(game.pause);
            }
            key = _getch();
            if (key == 0 || key == 224) {
                key = _getch();
                if (key == 72 && game.dir != Down) {
                    game.dir = Up;
                    game.dx = 0;
                    game.dy = -1;
                }
                else if (key == 80 && game.dir != Up) {
                    game.dir = Down;
                    game.dx = 0;
                    game.dy = 1;
                }
                else if (key == 75 && game.dir != Right) {
                    game.dir = Left;
                    game.dx = -1;
                    game.dy = 0;
                }
                else if (key == 77 && game.dir != Left) {
                    game.dir = Right;
                    game.dx = 1;
                    game.dy = 0;
                }
            }
        }
    }
}