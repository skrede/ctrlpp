#include "build_id.h"

#include <cstddef>
#include <cstdint>
#include <cstring>

extern "C" const std::uint8_t __start_note_gnu_build_id[];
extern "C" const std::uint8_t __stop_note_gnu_build_id[];

namespace
{

struct note_header
{
    std::uint32_t name_size;
    std::uint32_t descriptor_size;
    std::uint32_t type;
};

constexpr std::uint32_t kNoteTypeGnuBuildId = 3;
constexpr std::uint32_t kOwnerSize          = 4;
constexpr std::uint32_t kDescriptorSize     = 20;

// The note is a raw byte run delimited by two linker-provided symbols, so there
// is no object to bind a reference to and the layout the System V gABI gives it
// -- an Elf32_Nhdr, then the padded owner string, then the descriptor -- has to
// be walked through the addresses themselves. A null return is the modelled
// absence of a well-formed note.
const std::uint8_t *build_id_descriptor() noexcept
{
    const std::uint8_t *note = __start_note_gnu_build_id;
    const std::size_t   span = static_cast<std::size_t>(__stop_note_gnu_build_id - note);
    if(span < sizeof(note_header) + kOwnerSize + kDescriptorSize)
        return nullptr;

    note_header header{};
    std::memcpy(&header, note, sizeof(header));
    if(header.type != kNoteTypeGnuBuildId || header.name_size != kOwnerSize
       || header.descriptor_size != kDescriptorSize)
        return nullptr;

    if(std::memcmp(note + sizeof(header), "GNU", kOwnerSize) != 0)
        return nullptr;

    return note + sizeof(header) + kOwnerSize;
}

}

namespace ctrlpp
{

void render_build_id(char (&text)[kBuildIdTextLength]) noexcept
{
    const std::uint8_t *descriptor = build_id_descriptor();
    if(descriptor == nullptr)
    {
        std::memcpy(text, "unavailable", sizeof("unavailable"));
        return;
    }

    constexpr char digits[] = "0123456789abcdef";
    for(std::size_t i = 0; i < kDescriptorSize; ++i)
    {
        text[2 * i]     = digits[descriptor[i] >> 4];
        text[2 * i + 1] = digits[descriptor[i] & 0x0Fu];
    }
    text[2 * kDescriptorSize] = '\0';
}

}
