#ifndef TIMER_H_
#define TIMER_H_

#include <chrono>

// Using a namespace to prevent potential naming conflicts and to logically group
// associated functionalities.
namespace TimerUtils {

    /**
    * @class Timer
    * Represents a timer that allows the user to pause, resume, or return to the main menu.
    */
    class Timer {
    private:
        std::chrono::steady_clock::time_point startTime;
        std::chrono::steady_clock::duration pausedTime;
        bool isPaused;

    public:
        /**
        * Constructs a Timer object and initializes its properties.
        */
        Timer() {
            isPaused = false;
        }

        /**
        * Starts the timer.
        */
        void start() {
            startTime = std::chrono::steady_clock::now();
            pausedTime = std::chrono::steady_clock::duration::zero();
            isPaused = false;
        }

        /**
        * Pauses the timer.
        */
        void pause() {
            if (!isPaused) {
                pausedTime += std::chrono::steady_clock::now() - startTime;
                isPaused = true;
            }
        }

        /**
        * Resumes the timer.
        */
        void resume() {
            if (isPaused) {
                startTime = std::chrono::steady_clock::now();
                isPaused = false;
            }
        }

        /**
        * Returns the elapsed time in seconds.
        *
        * @return double The elapsed time in seconds.
        */
        double getElapsedTime() {
            std::chrono::steady_clock::duration elapsedTime;

            if (isPaused) {
                elapsedTime = pausedTime;
            }
            else {
                elapsedTime = std::chrono::steady_clock::now() - startTime + pausedTime;
            }

            return std::chrono::duration<double>(elapsedTime).count();
        }
    };
}

#endif  /* !TIMER_H_ */
