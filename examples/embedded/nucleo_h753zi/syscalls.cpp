// The system calls the C library's stdio and abort reach on this image.
// Descriptors 0, 1 and 2 are the USART3 console, a write-only character device,
// and nothing else is ever open; the image is the only process.

#include <span>
#include <array>
#include <cerrno>
#include <csignal>
#include <cstddef>
#include <cstdlib>
#include <unistd.h>

#include <sys/stat.h>

namespace ctrlpp {
void usart3_write(std::span<const char> text) noexcept;
}

namespace {

constexpr int kConsoleDescriptors = 3;
constexpr pid_t kImageProcess     = 1;

std::array<bool, kConsoleDescriptors> console_open{true, true, true};

bool is_console(int descriptor) noexcept
{
    return descriptor >= 0 && descriptor < kConsoleDescriptors && console_open[static_cast<std::size_t>(descriptor)];
}

int refuse(int code) noexcept
{
    errno = code;
    return -1;
}

}

extern "C" ssize_t _write(int descriptor, const void *data, std::size_t length)
{
    if(!is_console(descriptor))
        return refuse(EBADF);
    ctrlpp::usart3_write({static_cast<const char *>(data), length});
    return static_cast<ssize_t>(length);
}

extern "C" ssize_t _read(int descriptor, void *, std::size_t)
{
    return is_console(descriptor) ? 0 : refuse(EBADF);
}

extern "C" off_t _lseek(int descriptor, off_t, int)
{
    return refuse(is_console(descriptor) ? ESPIPE : EBADF);
}

extern "C" int _close(int descriptor)
{
    if(!is_console(descriptor))
        return refuse(EBADF);
    console_open[static_cast<std::size_t>(descriptor)] = false;
    return 0;
}

// The stream buffering setup reads the device type here and line-buffers a
// character device.
extern "C" int _fstat(int descriptor, struct stat *status)
{
    if(!is_console(descriptor))
        return refuse(EBADF);
    *status         = {};
    status->st_mode = S_IFCHR;
    return 0;
}

extern "C" int _isatty(int descriptor)
{
    if(is_console(descriptor))
        return 1;
    refuse(EBADF);
    return 0;
}

extern "C" pid_t _getpid()
{
    return kImageProcess;
}

// Every signal's action on this image is to stop it, as abort's SIGABRT does.
extern "C" int _kill(pid_t process, int signal)
{
    if(process != kImageProcess)
        return refuse(ESRCH);
    if(signal < 0 || signal >= NSIG)
        return refuse(EINVAL);
    if(signal == 0)
        return 0;
    _exit(EXIT_FAILURE);
}

extern "C" void _exit(int)
{
    for(;;)
    {
    }
}
