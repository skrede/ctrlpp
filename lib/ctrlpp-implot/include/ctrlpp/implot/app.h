#ifndef HPP_GUARD_CTRLPP_IMPLOT_APP_H
#define HPP_GUARD_CTRLPP_IMPLOT_APP_H

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "implot.h"

#include <GLFW/glfw3.h>

#include <concepts>

namespace ctrlpp::implot {

class App {
public:
    App(int width, int height, const char* title);
    ~App();

    App(const App&) = delete;
    App& operator=(const App&) = delete;
    App(App&&) = delete;
    App& operator=(App&&) = delete;

    template <std::invocable F>
    void run(F&& frame_fn)
    {
        while (!glfwWindowShouldClose(window_)) {
            glfwPollEvents();

            if (glfwGetWindowAttrib(window_, GLFW_ICONIFIED)) {
                glfwWaitEvents();
                continue;
            }

            ImGui_ImplOpenGL3_NewFrame();
            ImGui_ImplGlfw_NewFrame();
            ImGui::NewFrame();

            frame_fn();

            ImGui::Render();
            int display_w = 0;
            int display_h = 0;
            glfwGetFramebufferSize(window_, &display_w, &display_h);
            glViewport(0, 0, display_w, display_h);
            glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT);
            ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
            glfwSwapBuffers(window_);
        }
    }

    [[nodiscard]] bool should_close() const;

private:
    GLFWwindow* window_{nullptr};
};

}

#endif
